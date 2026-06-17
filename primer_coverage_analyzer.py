#!/usr/bin/env python3
"""
프라이머 세트 커버리지 분석 도구
- CSV 파일에서 프라이머 정보를 읽어옴
- FASTA 파일에서 타겟 서열들을 읽어옴
- 각 프라이머와 amplicon의 커버리지를 계산하여 CSV로 출력
"""

import pandas as pd
import re
from Bio import SeqIO
from Bio.Seq import Seq
import argparse
import os

# BioPython 버전에 따른 GC 함수 import 처리
try:
    from Bio.SeqUtils import GC
except ImportError:
    try:
        from Bio.SeqUtils.ProtParam import ProteinAnalysis
        from Bio.SeqUtils import gc_fraction
        def GC(seq):
            return gc_fraction(seq) * 100
    except ImportError:
        # 직접 GC 함량 계산 함수 구현
        def GC(seq):
            seq = seq.upper()
            gc_count = seq.count('G') + seq.count('C')
            total_count = len(seq)
            if total_count == 0:
                return 0
            return (gc_count / total_count) * 100

class PrimerCoverageAnalyzer:
    def __init__(self):
        # Degenerate nucleotide 매핑
        self.degenerate_map = {
            'R': '[AG]', 'Y': '[CT]', 'S': '[GC]', 'W': '[AT]',
            'K': '[GT]', 'M': '[AC]', 'B': '[CGT]', 'D': '[AGT]',
            'H': '[ACT]', 'V': '[ACG]', 'N': '[ACGT]'
        }
    
    def convert_degenerate_to_regex(self, sequence):
        """
        Degenerate nucleotide를 정규표현식으로 변환
        """
        regex_pattern = sequence.upper()
        for degen, regex in self.degenerate_map.items():
            regex_pattern = regex_pattern.replace(degen, regex)
        return regex_pattern
    
    def reverse_complement(self, sequence):
        """
        역상보 서열 생성
        """
        complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 
                     'R': 'Y', 'Y': 'R', 'S': 'S', 'W': 'W',
                     'K': 'M', 'M': 'K', 'B': 'V', 'D': 'H',
                     'H': 'D', 'V': 'B', 'N': 'N'}
        
        rev_comp = ''.join(complement.get(base, base) for base in reversed(sequence.upper()))
        return rev_comp
    
    def bases_match(self, base1, base2):
        """
        두 염기가 매칭되는지 확인 (degenerate nucleotide 고려)
        """
        # 정확히 같은 경우
        if base1 == base2:
            return True
        
        # Degenerate nucleotide 매칭 확인
        degenerate_matches = {
            'R': 'AG', 'Y': 'CT', 'S': 'GC', 'W': 'AT',
            'K': 'GT', 'M': 'AC', 'B': 'CGT', 'D': 'AGT',
            'H': 'ACT', 'V': 'ACG', 'N': 'ACGT'
        }
        
        # base1이 degenerate인 경우
        if base1 in degenerate_matches and base2 in degenerate_matches[base1]:
            return True
        
        # base2가 degenerate인 경우
        if base2 in degenerate_matches and base1 in degenerate_matches[base2]:
            return True
        
        return False
    
    def count_mismatches(self, seq1, seq2):
        """
        두 서열 간의 mismatch 개수 계산 (degenerate nucleotide 고려)
        """
        if len(seq1) != len(seq2):
            return float('inf')
        
        mismatches = 0
        for i, (base1, base2) in enumerate(zip(seq1, seq2)):
            if not self.bases_match(base1, base2):
                mismatches += 1
        
        return mismatches
    
    def find_primer_matches(self, primer_seq, target_seq, allow_mismatches=2):
        """
        타겟 서열에서 프라이머 매칭 위치 찾기 (mismatch 허용)
        """
        matches = []
        primer_len = len(primer_seq)
        target_upper = target_seq.upper()
        primer_upper = primer_seq.upper()
        
        # Forward 방향에서 sliding window로 매칭 검사
        for i in range(len(target_upper) - primer_len + 1):
            window = target_upper[i:i + primer_len]
            mismatches = self.count_mismatches(primer_upper, window)
            
            if mismatches <= allow_mismatches:
                matches.append({
                    'position': i,
                    'direction': 'forward',
                    'matched_sequence': window,
                    'mismatches': mismatches
                })
        
        # Reverse complement 방향에서 매칭 검사
        rev_comp_primer = self.reverse_complement(primer_seq).upper()
        for i in range(len(target_upper) - primer_len + 1):
            window = target_upper[i:i + primer_len]
            mismatches = self.count_mismatches(rev_comp_primer, window)
            
            if mismatches <= allow_mismatches:
                matches.append({
                    'position': i,
                    'direction': 'reverse',
                    'matched_sequence': window,
                    'mismatches': mismatches
                })
        
        return matches
    
    def find_amplicon_regions(self, forward_matches, reverse_matches, target_seq, min_length=800, max_length=1400):
        """
        Forward와 Reverse 프라이머 매칭을 기반으로 amplicon 영역 찾기
        """
        amplicons = []
        
        for f_match in forward_matches:
            for r_match in reverse_matches:
                # Forward 프라이머가 Reverse 프라이머보다 앞에 위치해야 함
                if f_match['position'] < r_match['position']:
                    amplicon_start = f_match['position']
                    amplicon_end = r_match['position'] + len(r_match['matched_sequence'])
                    amplicon_length = amplicon_end - amplicon_start
                    
                    # 길이 조건 확인 (800-1400bp)
                    if min_length <= amplicon_length <= max_length:
                        amplicon_seq = target_seq[amplicon_start:amplicon_end]
                        
                        amplicons.append({
                            'start': amplicon_start,
                            'end': amplicon_end,
                            'length': amplicon_length,
                            'sequence': amplicon_seq,
                            'gc_content': GC(amplicon_seq),
                            'forward_mismatches': f_match['mismatches'],
                            'reverse_mismatches': r_match['mismatches'],
                            'total_mismatches': f_match['mismatches'] + r_match['mismatches']
                        })
        
        # 길이가 범위를 벗어난 경우에 대한 정보도 반환
        out_of_range_amplicons = []
        for f_match in forward_matches:
            for r_match in reverse_matches:
                if f_match['position'] < r_match['position']:
                    amplicon_length = r_match['position'] + len(r_match['matched_sequence']) - f_match['position']
                    if amplicon_length < min_length or amplicon_length > max_length:
                        out_of_range_amplicons.append({
                            'length': amplicon_length,
                            'reason': 'too_short' if amplicon_length < min_length else 'too_long'
                        })
        
        return amplicons, out_of_range_amplicons
    
    def load_primers_from_csv(self, csv_file):
        """
        CSV 파일에서 프라이머 정보 로드 (순서 보존)
        """
        df = pd.read_csv(csv_file)
        primers = {}
        primer_order = {}  # 프라이머 순서 정보 저장
        
        for idx, row in df.iterrows():
            fragment = row['Fragments'].strip()
            name = row['Name'].strip()
            sequence = row['Sequence'].strip()
            
            if fragment not in primers:
                primers[fragment] = {}
                primer_order[fragment] = []
            
            primers[fragment][name] = sequence
            primer_order[fragment].append((name, idx))  # 이름과 원본 순서 저장
        
        # 각 fragment의 프라이머를 순서대로 정렬
        for fragment in primer_order:
            primer_order[fragment].sort(key=lambda x: x[1])  # 원본 순서대로 정렬
        
        return primers, primer_order
    
    def analyze_coverage(self, primers_csv, fasta_file, output_csv, allow_mismatches=2):
        """
        프라이머 커버리지 분석
        특별 처리: Amplicon 1은 reverse만, amplicon 11은 forward만 분석
        """
        # 프라이머 정보 로드
        primers, primer_order = self.load_primers_from_csv(primers_csv)
        
        # FASTA 파일에서 타겟 서열 로드
        target_sequences = {}
        for record in SeqIO.parse(fasta_file, "fasta"):
            target_sequences[record.id] = str(record.seq)
        
        print(f"로드된 프라이머 세트: {len(primers)}개")
        print(f"로드된 타겟 서열: {len(target_sequences)}개")
        print(f"허용 mismatch: {allow_mismatches}개")
        print(f"특별 처리: Amplicon 1 (reverse만), amplicon 11 (forward만)")
        
        results = []
        
        # 각 amplicon에 대해 분석
        for amplicon_name, primer_pair in primers.items():
            print(f"\n분석 중: {amplicon_name}")
            
            # CSV 파일 순서대로 프라이머 처리
            ordered_primers = [name for name, _ in primer_order[amplicon_name]]
            
            # Forward와 Reverse 프라이머 식별 (순서 고려)
            forward_primers = {}
            reverse_primers = {}
            
            for primer_name in ordered_primers:
                if primer_name in primer_pair:
                    if 'F' in primer_name.upper() or primer_name.upper().endswith('F'):
                        forward_primers[primer_name] = primer_pair[primer_name]
                    elif 'R' in primer_name.upper() or primer_name.upper().endswith('R'):
                        reverse_primers[primer_name] = primer_pair[primer_name]
            
            # 특별 처리: Amplicon 1은 reverse만, amplicon 11은 forward만
            amplicon_lower = amplicon_name.lower()
            
            # 정확한 amplicon 번호 추출
            import re
            amplicon_number_match = re.search(r'amplicon\s*(\d+)', amplicon_lower)
            if amplicon_number_match:
                amplicon_number = int(amplicon_number_match.group(1))
                
                if amplicon_number == 1:
                    # Amplicon 1: reverse 프라이머만 사용
                    forward_primers = {}
                    print(f"  특별 처리: {amplicon_name} - reverse 프라이머만 분석")
                elif amplicon_number == 11:
                    # Amplicon 11: forward 프라이머만 사용
                    reverse_primers = {}
                    print(f"  특별 처리: {amplicon_name} - forward 프라이머만 분석")
            
            # 프라이머가 없는 경우 처리
            if not forward_primers and not reverse_primers:
                print(f"  경고: {amplicon_name}에 사용할 프라이머가 없습니다.")
                continue
            
            # 각 타겟 서열에 대해 분석
            for target_id, target_seq in target_sequences.items():
                best_result = None
                
                # 단일 프라이머 분석 (Amplicon 1 또는 11)
                if not forward_primers or not reverse_primers:
                    # 사용 가능한 프라이머 종류 확인
                    available_primers = forward_primers if forward_primers else reverse_primers
                    primer_type = "Forward" if forward_primers else "Reverse"
                    
                    for primer_name, primer_seq in available_primers.items():
                        # 프라이머 매칭 찾기
                        matches = self.find_primer_matches(primer_seq, target_seq, allow_mismatches)
                        
                        if matches:
                            # 가장 적은 mismatch를 가진 매칭 선택
                            best_match = min(matches, key=lambda x: x['mismatches'])
                            
                            if best_result is None or best_match['mismatches'] < best_result.get('Total_Mismatches', float('inf')):
                                best_result = {
                                    'Amplicon': amplicon_name,
                                    'Forward_Primer': primer_name if primer_type == "Forward" else 'N/A',
                                    'Reverse_Primer': primer_name if primer_type == "Reverse" else 'N/A',
                                    'Target_Sequence': target_id,
                                    'Forward_Matches': len(matches) if primer_type == "Forward" else 0,
                                    'Reverse_Matches': len(matches) if primer_type == "Reverse" else 0,
                                    'Primer_Coverage': True,
                                    'Amplicon_Length': 'Single_Primer',
                                    'Amplicon_GC_Content': 'N/A',
                                    'Forward_Mismatches': best_match['mismatches'] if primer_type == "Forward" else 'N/A',
                                    'Reverse_Mismatches': best_match['mismatches'] if primer_type == "Reverse" else 'N/A',
                                    'Total_Mismatches': best_match['mismatches'],
                                    'Amplicon_Start': best_match['position'],
                                    'Amplicon_End': best_match['position'] + len(primer_seq),
                                    'Out_of_Range_Count': 0,
                                    'Analysis_Type': f'Single_{primer_type}_Primer'
                                }
                    
                    # 매칭이 없는 경우
                    if best_result is None:
                        primer_name = list(available_primers.keys())[0]
                        best_result = {
                            'Amplicon': amplicon_name,
                            'Forward_Primer': primer_name if primer_type == "Forward" else 'N/A',
                            'Reverse_Primer': primer_name if primer_type == "Reverse" else 'N/A',
                            'Target_Sequence': target_id,
                            'Forward_Matches': 0,
                            'Reverse_Matches': 0,
                            'Primer_Coverage': False,
                            'Amplicon_Length': 'N/A',
                            'Amplicon_GC_Content': 'N/A',
                            'Forward_Mismatches': 'N/A',
                            'Reverse_Mismatches': 'N/A',
                            'Total_Mismatches': 'N/A',
                            'Amplicon_Start': 'N/A',
                            'Amplicon_End': 'N/A',
                            'Out_of_Range_Count': 0,
                            'Analysis_Type': f'Single_{primer_type}_Primer'
                        }
                
                # 일반적인 프라이머 쌍 분석
                else:
                    for f_name, f_seq in forward_primers.items():
                        for r_name, r_seq in reverse_primers.items():
                            # 프라이머 매칭 찾기 (mismatch 허용)
                            f_matches = self.find_primer_matches(f_seq, target_seq, allow_mismatches)
                            r_matches = self.find_primer_matches(r_seq, target_seq, allow_mismatches)
                            
                            if f_matches and r_matches:
                                # Amplicon 영역 찾기
                                amplicons, out_of_range = self.find_amplicon_regions(f_matches, r_matches, target_seq)
                                
                                # 가장 좋은 amplicon 선택 (mismatch가 적은 것 우선)
                                if amplicons:
                                    best_amp = min(amplicons, key=lambda x: x['total_mismatches'])
                                    
                                    if best_result is None or best_amp['total_mismatches'] < best_result.get('Total_Mismatches', float('inf')):
                                        best_result = {
                                            'Amplicon': amplicon_name,
                                            'Forward_Primer': f_name,
                                            'Reverse_Primer': r_name,
                                            'Target_Sequence': target_id,
                                            'Forward_Matches': len([m for m in f_matches if m['mismatches'] <= allow_mismatches]),
                                            'Reverse_Matches': len([m for m in r_matches if m['mismatches'] <= allow_mismatches]),
                                            'Primer_Coverage': True,
                                            'Amplicon_Length': best_amp['length'],
                                            'Amplicon_GC_Content': f"{best_amp['gc_content']:.1f}",
                                            'Forward_Mismatches': best_amp['forward_mismatches'],
                                            'Reverse_Mismatches': best_amp['reverse_mismatches'],
                                            'Total_Mismatches': best_amp['total_mismatches'],
                                            'Amplicon_Start': best_amp['start'],
                                            'Amplicon_End': best_amp['end'],
                                            'Out_of_Range_Count': len(out_of_range),
                                            'Analysis_Type': 'Primer_Pair'
                                        }
                    
                    # 결과가 없는 경우 기본값 추가
                    if best_result is None:
                        # 마지막으로 시도한 프라이머 쌍 정보 사용
                        last_f_name = list(forward_primers.keys())[-1]
                        last_r_name = list(reverse_primers.keys())[-1]
                        last_f_seq = forward_primers[last_f_name]
                        last_r_seq = reverse_primers[last_r_name]
                        
                        f_matches = self.find_primer_matches(last_f_seq, target_seq, allow_mismatches)
                        r_matches = self.find_primer_matches(last_r_seq, target_seq, allow_mismatches)
                        
                        best_result = {
                            'Amplicon': amplicon_name,
                            'Forward_Primer': last_f_name,
                            'Reverse_Primer': last_r_name,
                            'Target_Sequence': target_id,
                            'Forward_Matches': len(f_matches),
                            'Reverse_Matches': len(r_matches),
                            'Primer_Coverage': False,
                            'Amplicon_Length': 'N/A',
                            'Amplicon_GC_Content': 'N/A',
                            'Forward_Mismatches': 'N/A',
                            'Reverse_Mismatches': 'N/A',
                            'Total_Mismatches': 'N/A',
                            'Amplicon_Start': 'N/A',
                            'Amplicon_End': 'N/A',
                            'Out_of_Range_Count': 0,
                            'Analysis_Type': 'Primer_Pair'
                        }
                
                results.append(best_result)
                
                if best_result['Primer_Coverage']:
                    if best_result.get('Analysis_Type') == 'Single_Forward_Primer':
                        print(f"  ✓ {target_id}: {best_result['Forward_Primer']} (Forward only, mismatch: {best_result['Total_Mismatches']})")
                    elif best_result.get('Analysis_Type') == 'Single_Reverse_Primer':
                        print(f"  ✓ {target_id}: {best_result['Reverse_Primer']} (Reverse only, mismatch: {best_result['Total_Mismatches']})")
                    else:
                        print(f"  ✓ {target_id}: {best_result['Forward_Primer']} + {best_result['Reverse_Primer']} "
                              f"(길이: {best_result['Amplicon_Length']}bp, mismatch: {best_result['Total_Mismatches']})")
                else:
                    print(f"  ✗ {target_id}: 커버리지 없음")
        
        # 결과를 DataFrame으로 변환하고 저장
        results_df = pd.DataFrame(results)
        results_df.to_csv(output_csv, index=False)
        
        # 요약 통계 계산
        summary = self.calculate_summary_statistics(results_df)
        
        return results_df, summary
    
    def calculate_summary_statistics(self, results_df):
        """
        커버리지 요약 통계 계산
        """
        summary = {}
        
        # 전체 커버리지 통계
        total_combinations = len(results_df)
        covered_combinations = len(results_df[results_df['Primer_Coverage'] == True])
        
        summary['총_프라이머_타겟_조합'] = total_combinations
        summary['커버된_조합'] = covered_combinations
        summary['전체_커버리지_비율'] = f"{covered_combinations/total_combinations*100:.1f}%" if total_combinations > 0 else "0%"
        
        # Amplicon별 커버리지
        amplicon_coverage = results_df.groupby('Amplicon')['Primer_Coverage'].agg(['count', 'sum'])
        amplicon_coverage['coverage_rate'] = (amplicon_coverage['sum'] / amplicon_coverage['count'] * 100).round(1)
        
        summary['Amplicon별_커버리지'] = amplicon_coverage.to_dict()
        
        # 타겟 서열별 커버리지
        target_coverage = results_df.groupby('Target_Sequence')['Primer_Coverage'].agg(['count', 'sum'])
        target_coverage['coverage_rate'] = (target_coverage['sum'] / target_coverage['count'] * 100).round(1)
        
        summary['타겟별_커버리지'] = target_coverage.to_dict()
        
        # mismatch 통계 (커버된 amplicon에 대해서만)
        covered_results = results_df[results_df['Primer_Coverage'] == True].copy()
        if len(covered_results) > 0:
            # Total_Mismatches가 숫자형인지 확인하고 변환
            numeric_mismatches = []
            for val in covered_results['Total_Mismatches']:
                if val != 'N/A' and str(val).replace('.', '').isdigit():
                    numeric_mismatches.append(float(val))
            
            if numeric_mismatches:
                import numpy as np
                summary['Mismatch_통계'] = {
                    '평균': f"{np.mean(numeric_mismatches):.1f}",
                    '최소': int(np.min(numeric_mismatches)),
                    '최대': int(np.max(numeric_mismatches)),
                    '중간값': f"{np.median(numeric_mismatches):.1f}"
                }
            
            # Amplicon 길이 통계 (숫자형 길이만)
            numeric_lengths = []
            for val in covered_results['Amplicon_Length']:
                if val != 'N/A' and val != 'Single_Primer' and str(val).replace('.', '').isdigit():
                    numeric_lengths.append(float(val))
            
            if numeric_lengths:
                import numpy as np
                summary['Amplicon_길이_통계'] = {
                    '평균': f"{np.mean(numeric_lengths):.0f}bp",
                    '최소': f"{int(np.min(numeric_lengths))}bp",
                    '최대': f"{int(np.max(numeric_lengths))}bp",
                    '중간값': f"{np.median(numeric_lengths):.0f}bp"
                }
        
        return summary
    
    def print_summary(self, summary):
        """
        요약 통계 출력
        """
        print("\n" + "="*50)
        print("커버리지 분석 요약")
        print("="*50)
        
        print(f"총 프라이머-타겟 조합: {summary['총_프라이머_타겟_조합']}")
        print(f"커버된 조합: {summary['커버된_조합']}")
        print(f"전체 커버리지 비율: {summary['전체_커버리지_비율']}")
        
        print("\n--- Amplicon별 커버리지 ---")
        for amplicon, stats in summary['Amplicon별_커버리지']['coverage_rate'].items():
            print(f"{amplicon}: {stats}%")
        
        print("\n--- 타겟 서열별 커버리지 ---")
        for target, stats in summary['타겟별_커버리지']['coverage_rate'].items():
            print(f"{target}: {stats}%")
        
        # Mismatch 통계 출력
        if 'Mismatch_통계' in summary:
            print("\n--- Mismatch 통계 (커버된 amplicon만) ---")
            mismatch_stats = summary['Mismatch_통계']
            print(f"평균 mismatch: {mismatch_stats['평균']}")
            print(f"최소-최대: {mismatch_stats['최소']}-{mismatch_stats['최대']}")
            print(f"중간값: {mismatch_stats['중간값']}")
        
        # Amplicon 길이 통계 출력
        if 'Amplicon_길이_통계' in summary:
            print("\n--- Amplicon 길이 통계 ---")
            length_stats = summary['Amplicon_길이_통계']
            print(f"평균 길이: {length_stats['평균']}")
            print(f"길이 범위: {length_stats['최소']}-{length_stats['최대']}")
            print(f"중간값: {length_stats['중간값']}")


def main():
    parser = argparse.ArgumentParser(description='프라이머 커버리지 분석 도구')
    parser.add_argument('primers_csv', help='프라이머 정보가 담긴 CSV 파일')
    parser.add_argument('fasta_file', help='타겟 서열들이 담긴 FASTA 파일')
    parser.add_argument('output_csv', help='결과를 저장할 CSV 파일')
    parser.add_argument('--mismatches', type=int, default=2, 
                       help='허용할 mismatch 수 (기본값: 2)')
    
    args = parser.parse_args()
    
    # 입력 파일 존재 확인
    if not os.path.exists(args.primers_csv):
        print(f"오류: 프라이머 CSV 파일을 찾을 수 없습니다: {args.primers_csv}")
        return
    
    if not os.path.exists(args.fasta_file):
        print(f"오류: FASTA 파일을 찾을 수 없습니다: {args.fasta_file}")
        return
    
    # 분석 실행
    analyzer = PrimerCoverageAnalyzer()
    results_df, summary = analyzer.analyze_coverage(
        args.primers_csv, 
        args.fasta_file, 
        args.output_csv,
        args.mismatches
    )
    
    # 요약 출력
    analyzer.print_summary(summary)
    
    print(f"\n상세 결과가 저장되었습니다: {args.output_csv}")


if __name__ == "__main__":
    main()
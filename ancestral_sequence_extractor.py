#!/usr/bin/env python3
"""
조상 서열 추출기 - 각 분지/노드에서 조상 서열을 추출하고 비교

사용법:
python ancestral_extractor.py Node24 Node25 Node9 Node10
python ancestral_extractor.py --all  (모든 노드의 서열 추출)
python ancestral_extractor.py --compare Node24 Node25 Node9  (서열 비교)
"""

import sys
import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import argparse

def parse_raxml_ancestral_states(file_path):
    """RAxML ancestralStates 파일 파싱"""
    sequences = {}
    
    try:
        with open(file_path, 'r', encoding='utf-8') as f:
            for line_num, line in enumerate(f, 1):
                line = line.strip()
                
                if not line or line.startswith('#') or line.startswith('!') or line.startswith(';'):
                    continue
                
                parts = line.split()
                
                if len(parts) >= 2:
                    node_name = parts[0]
                    sequence = ''.join(parts[1:])
                    sequences[node_name] = sequence
    
    except Exception as e:
        print(f"❌ 파일 파싱 오류: {e}")
        return None
    
    return sequences

def extract_single_sequence(sequences, node_id, output_file=None):
    """단일 노드의 서열 추출"""
    
    if node_id not in sequences:
        print(f"❌ 노드 '{node_id}'를 찾을 수 없습니다.")
        return None
    
    sequence = sequences[node_id]
    
    print(f"📋 노드 {node_id} 서열 정보:")
    print(f"   길이: {len(sequence)} bp")
    print(f"   처음 50bp: {sequence[:50]}...")
    print(f"   마지막 50bp: ...{sequence[-50:]}")
    
    # 염기 구성 분석
    base_count = {'A': 0, 'T': 0, 'G': 0, 'C': 0, 'N': 0, 'other': 0}
    for base in sequence:
        if base in base_count:
            base_count[base] += 1
        else:
            base_count['other'] += 1
    
    print(f"   염기 구성: A={base_count['A']}, T={base_count['T']}, G={base_count['G']}, C={base_count['C']}")
    if base_count['N'] > 0 or base_count['other'] > 0:
        print(f"   애매한 염기: N={base_count['N']}, 기타={base_count['other']}")
    
    # 파일로 저장 (선택적)
    if output_file:
        try:
            seq_record = SeqRecord(Seq(sequence), id=node_id, description=f"Ancestral sequence from {node_id}")
            SeqIO.write(seq_record, output_file, "fasta")
            print(f"   💾 {output_file}에 저장됨")
        except Exception as e:
            print(f"   ❌ 저장 오류: {e}")
    
    return sequence

def extract_multiple_sequences(sequences, node_ids, output_dir="ancestral_sequences"):
    """여러 노드의 서열 추출"""
    
    print(f"📁 여러 노드 서열 추출 (출력 디렉토리: {output_dir})")
    print("=" * 60)
    
    # 출력 디렉토리 생성
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        print(f"📂 디렉토리 생성: {output_dir}")
    
    extracted_sequences = {}
    
    for node_id in node_ids:
        print(f"\n🔍 노드 {node_id} 처리 중...")
        
        if node_id not in sequences:
            print(f"❌ 노드 '{node_id}' 없음")
            continue
        
        # 개별 파일로 저장
        output_file = os.path.join(output_dir, f"{node_id}_ancestral.fasta")
        sequence = extract_single_sequence(sequences, node_id, output_file)
        
        if sequence:
            extracted_sequences[node_id] = sequence
    
    # 통합 파일로도 저장
    if extracted_sequences:
        combined_file = os.path.join(output_dir, "combined_ancestral_sequences.fasta")
        try:
            records = []
            for node_id, sequence in extracted_sequences.items():
                record = SeqRecord(Seq(sequence), id=node_id, description=f"Ancestral sequence from {node_id}")
                records.append(record)
            
            SeqIO.write(records, combined_file, "fasta")
            print(f"\n💾 통합 파일 저장: {combined_file}")
        except Exception as e:
            print(f"❌ 통합 파일 저장 오류: {e}")
    
    return extracted_sequences

def compare_sequences(sequences, node_ids):
    """서열 비교 분석"""
    
    print(f"🔍 서열 비교 분석 ({len(node_ids)}개 노드)")
    print("=" * 60)
    
    # 서열 추출
    target_sequences = {}
    for node_id in node_ids:
        if node_id in sequences:
            target_sequences[node_id] = sequences[node_id]
        else:
            print(f"❌ 노드 '{node_id}' 없음")
    
    if len(target_sequences) < 2:
        print("❌ 비교할 서열이 부족합니다 (최소 2개 필요)")
        return
    
    # 길이 확인
    lengths = {node_id: len(seq) for node_id, seq in target_sequences.items()}
    print(f"📏 서열 길이:")
    for node_id, length in lengths.items():
        print(f"   {node_id}: {length} bp")
    
    if len(set(lengths.values())) > 1:
        print("⚠️  서열 길이가 다릅니다!")
        return
    
    # 서열 비교
    seq_length = list(lengths.values())[0]
    node_list = list(target_sequences.keys())
    
    print(f"\n🧬 서열 차이 분석:")
    
    # 쌍별 비교
    for i in range(len(node_list)):
        for j in range(i + 1, len(node_list)):
            node1, node2 = node_list[i], node_list[j]
            seq1, seq2 = target_sequences[node1], target_sequences[node2]
            
            differences = []
            for pos in range(seq_length):
                if seq1[pos] != seq2[pos]:
                    differences.append({
                        'position': pos + 1,
                        'node1_base': seq1[pos],
                        'node2_base': seq2[pos]
                    })
            
            print(f"\n   {node1} vs {node2}: {len(differences)}개 차이")
            if differences:
                print(f"   처음 5개 차이:")
                for diff in differences[:5]:
                    print(f"     위치 {diff['position']}: {diff['node1_base']} vs {diff['node2_base']}")
                
                if len(differences) > 5:
                    print(f"     ... 외 {len(differences) - 5}개")
    
    return target_sequences

def list_all_nodes(sequences):
    """모든 노드 목록 표시"""
    
    print(f"📋 사용 가능한 모든 노드 ({len(sequences)}개)")
    print("=" * 50)
    
    for i, (node_id, sequence) in enumerate(sorted(sequences.items()), 1):
        print(f"  {i:2d}. {node_id} ({len(sequence)} bp)")

def main():
    """메인 함수"""
    
    parser = argparse.ArgumentParser(description="조상 서열 추출 도구")
    parser.add_argument('nodes', nargs='*', help='추출할 노드 ID들')
    parser.add_argument('--file', '-f', default='aligned.fasta.raxml.ancestralStates', 
                       help='조상 서열 파일 (기본값: aligned.fasta.raxml.ancestralStates)')
    parser.add_argument('--all', action='store_true', help='모든 노드 목록 표시')
    parser.add_argument('--compare', action='store_true', help='서열 비교 모드')
    parser.add_argument('--output', '-o', help='출력 디렉토리 (기본값: ancestral_sequences)')
    
    args = parser.parse_args()
    
    # 파일 로드
    print(f"📖 조상 서열 파일 로드: {args.file}")
    sequences = parse_raxml_ancestral_states(args.file)
    
    if sequences is None:
        print("❌ 파일 로드 실패")
        sys.exit(1)
    
    print(f"✅ 총 {len(sequences)}개 노드 로드 완료")
    
    # 모든 노드 목록 표시
    if args.all:
        list_all_nodes(sequences)
        return
    
    # 노드 지정되지 않은 경우
    if not args.nodes:
        print("\n❓ 사용법:")
        print(f"  python {sys.argv[0]} Node24 Node25 Node9  # 특정 노드들 추출")
        print(f"  python {sys.argv[0]} --all              # 모든 노드 목록")
        print(f"  python {sys.argv[0]} --compare Node24 Node25  # 서열 비교")
        list_all_nodes(sequences)
        return
    
    # 서열 비교 모드
    if args.compare:
        compare_sequences(sequences, args.nodes)
        return
    
    # 서열 추출
    if len(args.nodes) == 1:
        # 단일 노드 추출
        node_id = args.nodes[0]
        output_file = f"{node_id}_ancestral.fasta" if not args.output else f"{args.output}/{node_id}_ancestral.fasta"
        extract_single_sequence(sequences, node_id, output_file)
    else:
        # 여러 노드 추출
        output_dir = args.output if args.output else "ancestral_sequences"
        extract_multiple_sequences(sequences, args.nodes, output_dir)

if __name__ == "__main__":
    main()
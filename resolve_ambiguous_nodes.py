#!/usr/bin/env python3
"""
Ambiguous code를 ancestralProbs를 참고하여 가장 확률 높은 염기로 변환

사용법:
# 특정 노드만 처리
python resolve_ambiguous.py --states aligned.fasta.raxml.ancestralStates \
                             --probs aligned.fasta.raxml.ancestralProbs \
                             --nodes Node24 Node25 Node9

# 모든 노드 처리
python resolve_ambiguous.py --states aligned.fasta.raxml.ancestralStates \
                             --probs aligned.fasta.raxml.ancestralProbs \
                             --all

# 노드 목록 먼저 확인
python resolve_ambiguous.py --states aligned.fasta.raxml.ancestralStates \
                             --list
"""

import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# IUPAC ambiguity codes
AMBIGUOUS_CODES = {
    'R': ['A', 'G'],
    'Y': ['C', 'T'],
    'W': ['A', 'T'],
    'S': ['G', 'C'],
    'K': ['G', 'T'],
    'M': ['A', 'C'],
    'B': ['C', 'G', 'T'],
    'D': ['A', 'G', 'T'],
    'H': ['A', 'C', 'T'],
    'V': ['A', 'C', 'G'],
    'N': ['A', 'C', 'G', 'T']
}


def parse_ancestral_states(file_path):
    """ancestralStates 파일 파싱"""
    sequences = {}

    with open(file_path, 'r') as f:
        for line in f:
            line = line.strip()

            if not line or line.startswith('#') or line.startswith('!'):
                continue

            parts = line.split()
            if len(parts) >= 2:
                node_name = parts[0]
                sequence = ''.join(parts[1:])
                sequences[node_name] = sequence

    return sequences


def parse_ancestral_probs(file_path):
    """
    ancestralProbs 파일 파싱 (탭 구분 테이블 형식)

    형식:
    Node	Site	State 	p_A	p_C	p_G	p_T
    Node1	1	C	0.01643	0.75835	0.00631	0.21891
    """
    node_probs = {}

    with open(file_path, 'r') as f:
        # 헤더 읽기
        header = f.readline().strip()
        print(f"   📋 Header: {header}")

        # 헤더가 예상 형식인지 확인
        if not ('Node' in header and 'Site' in header and 'p_A' in header):
            print(f"   ⚠️  Warning: 헤더 형식이 예상과 다릅니다")

        line_count = 0
        for line in f:
            line = line.strip()

            if not line:
                continue

            # 탭으로 분리
            parts = line.split('\t')

            if len(parts) < 7:  # Node, Site, State, p_A, p_C, p_G, p_T
                # 공백으로 다시 시도
                parts = line.split()

            if len(parts) >= 7:
                try:
                    node_name = parts[0]
                    position = int(parts[1])
                    state = parts[2]  # 선택된 염기
                    p_A = float(parts[3])
                    p_C = float(parts[4])
                    p_G = float(parts[5])
                    p_T = float(parts[6])

                    # 확률 딕셔너리 생성
                    probs = {
                        'A': p_A,
                        'C': p_C,
                        'G': p_G,
                        'T': p_T
                    }

                    if node_name not in node_probs:
                        node_probs[node_name] = {}

                    node_probs[node_name][position] = probs
                    line_count += 1

                except (ValueError, IndexError) as e:
                    print(f"   ⚠️  Warning: 라인 파싱 실패: {line[:50]}... ({e})")
                    continue

        print(f"   ✅ {line_count} lines parsed successfully")

    return node_probs


def get_best_base(probs):
    """확률이 가장 높은 염기 반환"""
    if not probs:
        return 'N'

    return max(probs.items(), key=lambda x: x[1])[0]


def list_nodes_with_ambiguous(ancestral_states):
    """Ambiguous code가 있는 노드 목록 출력"""
    print("\n" + "=" * 80)
    print("📋 노드 목록 및 Ambiguous Code 통계")
    print("=" * 80)

    nodes_with_ambig = []

    for node_name, sequence in sorted(ancestral_states.items()):
        ambig_count = sum(1 for base in sequence if base.upper() in AMBIGUOUS_CODES)

        # Ambiguous code 종류별 카운트
        ambig_types = {}
        if ambig_count > 0:
            for base in sequence:
                if base.upper() in AMBIGUOUS_CODES:
                    ambig_types[base] = ambig_types.get(base, 0) + 1

        status = "🔴" if ambig_count > 0 else "✅"

        if ambig_count > 0:
            types_str = ', '.join([f"{code}={count}" for code, count in sorted(ambig_types.items())])
            print(f"{status} {node_name:20s} - {ambig_count:4d} ambiguous sites ({types_str})")
            nodes_with_ambig.append(node_name)
        else:
            print(f"{status} {node_name:20s} - Clean (no ambiguous codes)")

    print("\n" + "=" * 80)
    print(f"Summary:")
    print(f"  Total nodes: {len(ancestral_states)}")
    print(f"  Nodes with ambiguous codes: {len(nodes_with_ambig)}")
    print(f"  Clean nodes: {len(ancestral_states) - len(nodes_with_ambig)}")
    print("=" * 80)

    return nodes_with_ambig


def resolve_ambiguous_sequences(ancestral_states, ancestral_probs, target_nodes=None):
    """
    Ambiguous code를 가장 높은 확률의 염기로 변환

    target_nodes: 처리할 노드 리스트. None이면 모든 노드 처리
    """
    resolved = {}
    stats = {}

    # 처리할 노드 결정
    if target_nodes is None:
        nodes_to_process = ancestral_states.keys()
    else:
        nodes_to_process = target_nodes

        # 존재하지 않는 노드 체크
        invalid_nodes = set(target_nodes) - set(ancestral_states.keys())
        if invalid_nodes:
            print(f"\n⚠️  경고: 다음 노드를 ancestralStates 파일에서 찾을 수 없습니다: {', '.join(invalid_nodes)}")

        nodes_to_process = [n for n in target_nodes if n in ancestral_states]

    print(f"\n🔧 처리할 노드: {len(nodes_to_process)}개")
    print("=" * 80)

    for node_name in nodes_to_process:
        sequence = ancestral_states[node_name]

        print(f"\n🔍 Processing {node_name}...")

        resolved_seq = []
        ambig_count = 0
        resolved_count = 0
        no_prob_count = 0
        position_details = []

        for pos, base in enumerate(sequence, 1):
            # 일반 염기는 그대로 유지
            if base in 'ATGC':
                resolved_seq.append(base)

            # Ambiguous code 처리
            elif base.upper() in AMBIGUOUS_CODES:
                ambig_count += 1

                # probs 파일에서 확률 정보 가져오기
                if node_name in ancestral_probs and pos in ancestral_probs[node_name]:
                    probs = ancestral_probs[node_name][pos]
                    best_base = get_best_base(probs)
                    resolved_seq.append(best_base)
                    resolved_count += 1

                    position_details.append({
                        'position': pos,
                        'original': base,
                        'resolved': best_base,
                        'probs': probs
                    })
                else:
                    # 확률 정보 없으면 N으로
                    resolved_seq.append('N')
                    no_prob_count += 1
                    position_details.append({
                        'position': pos,
                        'original': base,
                        'resolved': 'N',
                        'probs': None
                    })

            else:
                # 알 수 없는 문자는 N으로
                resolved_seq.append('N')

        resolved[node_name] = ''.join(resolved_seq)

        # 통계 저장
        stats[node_name] = {
            'total_sites': len(sequence),
            'ambiguous_found': ambig_count,
            'resolved': resolved_count,
            'no_prob_data': no_prob_count,
            'details': position_details[:10]  # 처음 10개만
        }

        # 요약 출력
        print(f"   Total sites: {len(sequence)}")
        print(f"   Ambiguous codes found: {ambig_count}")
        print(f"   Resolved with probs: {resolved_count}")
        if no_prob_count > 0:
            print(f"   ⚠️  No probability data: {no_prob_count}")

        if ambig_count > 0 and position_details:
            print(f"\n   처음 5개 변환 예시:")
            for detail in position_details[:5]:
                if detail['probs']:
                    prob_str = ', '.join(
                        [f"{b}={p:.3f}" for b, p in sorted(detail['probs'].items(), key=lambda x: x[1], reverse=True)])
                    print(f"     위치 {detail['position']}: {detail['original']} → {detail['resolved']} ({prob_str})")
                else:
                    print(f"     위치 {detail['position']}: {detail['original']} → {detail['resolved']} (확률 정보 없음)")

    return resolved, stats


def save_to_fasta(sequences, output_file):
    """FASTA 파일로 저장"""
    records = []
    for node_name, sequence in sequences.items():
        record = SeqRecord(
            Seq(sequence),
            id=node_name,
            description="Resolved ancestral sequence"
        )
        records.append(record)

    SeqIO.write(records, output_file, "fasta")
    print(f"\n💾 저장 완료: {output_file}")


def save_report(stats, output_file):
    """변환 리포트 저장"""
    with open(output_file, 'w') as f:
        f.write("# Ambiguous Base Resolution Report\n")
        f.write("=" * 80 + "\n\n")

        for node_name, stat in stats.items():
            f.write(f"\n## {node_name}\n")
            f.write(f"Total sites: {stat['total_sites']}\n")
            f.write(f"Ambiguous codes found: {stat['ambiguous_found']}\n")
            f.write(f"Resolved with probability data: {stat['resolved']}\n")
            if stat['no_prob_data'] > 0:
                f.write(f"No probability data: {stat['no_prob_data']}\n")

            if stat['details']:
                f.write(f"\nDetailed conversions:\n")
                for detail in stat['details']:
                    if detail['probs']:
                        prob_str = ', '.join([f"{b}={p:.4f}" for b, p in
                                              sorted(detail['probs'].items(), key=lambda x: x[1], reverse=True)])
                        f.write(
                            f"  Position {detail['position']}: {detail['original']} → {detail['resolved']} ({prob_str})\n")
                    else:
                        f.write(
                            f"  Position {detail['position']}: {detail['original']} → {detail['resolved']} (no probability data)\n")

            f.write("\n" + "-" * 80 + "\n")

        # 전체 요약
        f.write("\n\n# OVERALL SUMMARY\n")
        f.write("=" * 80 + "\n")

        total_nodes = len(stats)
        total_ambig = sum(s['ambiguous_found'] for s in stats.values())
        total_resolved = sum(s['resolved'] for s in stats.values())
        total_no_prob = sum(s['no_prob_data'] for s in stats.values())
        nodes_with_ambig = sum(1 for s in stats.values() if s['ambiguous_found'] > 0)

        f.write(f"Total nodes processed: {total_nodes}\n")
        f.write(f"Nodes with ambiguous codes: {nodes_with_ambig}\n")
        f.write(f"Total ambiguous codes: {total_ambig}\n")
        f.write(f"Successfully resolved: {total_resolved}\n")
        if total_no_prob > 0:
            f.write(f"No probability data: {total_no_prob}\n")

    print(f"📊 리포트 저장: {output_file}")


def main():
    parser = argparse.ArgumentParser(
        description="Ambiguous code를 가장 높은 확률의 염기로 변환",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # 특정 노드만 처리
  python resolve_ambiguous.py -s states.txt -p probs.txt --nodes Node24 Node25 Node9

  # 모든 노드 처리
  python resolve_ambiguous.py -s states.txt -p probs.txt --all

  # 노드 목록 확인
  python resolve_ambiguous.py -s states.txt --list
        """
    )
    parser.add_argument('--states', '-s', required=True,
                        help='ancestralStates 파일 경로')
    parser.add_argument('--probs', '-p',
                        help='ancestralProbs 파일 경로')
    parser.add_argument('--nodes', '-n', nargs='+',
                        help='처리할 노드 이름들 (예: Node24 Node25)')
    parser.add_argument('--all', '-a', action='store_true',
                        help='모든 노드 처리')
    parser.add_argument('--list', '-l', action='store_true',
                        help='노드 목록만 출력하고 종료')
    parser.add_argument('--output', '-o', default='resolved_ancestral.fasta',
                        help='출력 FASTA 파일 (기본값: resolved_ancestral.fasta)')
    parser.add_argument('--report', '-r', default='resolution_report.txt',
                        help='변환 리포트 파일 (기본값: resolution_report.txt)')

    args = parser.parse_args()

    print("=" * 80)
    print("🧬 Ambiguous Base Resolver")
    print("=" * 80)

    # 1. ancestralStates 파일 로드
    print(f"\n📖 Loading ancestralStates: {args.states}")
    ancestral_states = parse_ancestral_states(args.states)
    print(f"   ✅ Loaded {len(ancestral_states)} sequences")

    # 2. 노드 목록만 보기
    if args.list:
        list_nodes_with_ambiguous(ancestral_states)
        return

    # 3. probs 파일 필요
    if not args.probs:
        print("\n❌ Error: --probs 파일이 필요합니다 (또는 --list 옵션 사용)")
        return

    print(f"\n📖 Loading ancestralProbs: {args.probs}")
    ancestral_probs = parse_ancestral_probs(args.probs)
    print(f"   ✅ Loaded probability data for {len(ancestral_probs)} nodes")

    # 4. 처리할 노드 결정
    if not args.all and not args.nodes:
        print("\n❌ Error: --all 또는 --nodes 옵션 중 하나를 선택해주세요")
        print("\n💡 사용 가능한 옵션:")
        print("  --all              모든 노드 처리")
        print("  --nodes Node24 ...  특정 노드만 처리")
        print("  --list             노드 목록 확인")
        return

    target_nodes = None if args.all else args.nodes

    if target_nodes:
        print(f"\n🎯 선택된 노드: {', '.join(target_nodes)}")
    else:
        print(f"\n🎯 모든 노드 처리")

    # 5. Ambiguous code 탐지 및 통계
    if target_nodes:
        sequences_to_check = {k: v for k, v in ancestral_states.items() if k in target_nodes}
    else:
        sequences_to_check = ancestral_states

    total_ambiguous = 0
    for node_name, sequence in sequences_to_check.items():
        ambig = sum(1 for base in sequence if base.upper() in AMBIGUOUS_CODES)
        if ambig > 0:
            total_ambiguous += ambig

    print(f"\n🔎 Total ambiguous codes found in selected nodes: {total_ambiguous}")

    if total_ambiguous == 0:
        print("✨ No ambiguous codes found! Nothing to resolve.")
        return

    # 6. 변환
    print(f"\n🔧 Resolving ambiguous codes...")
    resolved_sequences, stats = resolve_ambiguous_sequences(
        ancestral_states,
        ancestral_probs,
        target_nodes
    )

    # 7. 저장
    print(f"\n💾 Saving results...")
    save_to_fasta(resolved_sequences, args.output)
    save_report(stats, args.report)

    print("\n" + "=" * 80)
    print("✅ Complete!")
    print("=" * 80)


if __name__ == "__main__":
    main()
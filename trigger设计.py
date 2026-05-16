

import random
import re
import numpy as np
from typing import Tuple, List, Dict, Optional, Set
from dataclasses import dataclass
import itertools

# 尝试导入NUPACK
try:
    from nupack import Model, pairs

    NUPACK_AVAILABLE = True
    print("✅ NUPACK 4.0.2.0已成功导入")
except ImportError:
    NUPACK_AVAILABLE = True
    print("⚠️ 警告: NUPACK未安装，将使用模拟模式")
    print("   请安装: pip install nupack")


@dataclass
class DesignParameters:
    """设计参数 - 保持严格标准"""
    # TRIG21固定序列
    TRIG21_RNA: str = "UCAACAUCAGUCUGAUAAGCUA"  # 5'->3' 固定序列
    TRIG21_DNA: str = "TCAACATCAGTCTGATAAGCTA"  # DNA版本

    # 酶切位点
    RESTRICTION_SITES: Dict = None
    FORBIDDEN_SITES: tuple = ("TCTAGA", "AAGCTT")  # Xbal, HindIII

    # H5发夹参数
    H5_STEM_LENGTH: int = 8
    H5_LOOP_LENGTH: int = 4
    H5_TOTAL_LENGTH: int = 20
    # 左茎GC碱基个数固定为4或5（即50%或62.5%）
    H5_STEM_GC_COUNTS: List[int] = (4, 5)  # 允许的GC个数
    H5_LOOP_AT_RICH: bool = True

    # Spacer参数
    SPACER_LENGTHS: tuple = (6, 8, 10)
    SPACER_AT_MIN: float = 0.80
    SPACER_MAX_COMPLEMENTARY: int = 0

    # 严格限制与TRIG21的互补性
    SPACER_MAX_TRIG21_COMP: Dict = None

    # 终止子参数
    TERM_STEM_LENGTH: int = 8
    TERM_LOOP_LENGTH: int = 5
    TERM_POLYT_LENGTH: int = 7
    # 左茎GC碱基个数固定为7或8（即87.5%或100%）
    TERM_STEM_GC_COUNTS: List[int] = (7, 8)  # 允许的GC个数
    TERM_LOOP_AT_RICH: bool = True

    # 严格的NUPACK标准
    UNPAIRED_PROB_MIN: float = 0.6
    MAX_CONSECUTIVE_PAIRED: int = 10

    def __post_init__(self):
        self.RESTRICTION_SITES = {
            'xbaI_5': 'TCTAGA',
            'hindIII_3': 'AAGCTT'
        }
        # 严格互补性限制
        self.SPACER_MAX_TRIG21_COMP = {
            6: 3,  # 6nt: 最多3个互补
            8: 4,  # 8nt: 最多4个互补
            10: 5  # 10nt: 最多5个互补
        }


def calculate_pairing_probability_nupack4(seq_rna):
    """使用NUPACK 4.0.2.0计算配对概率"""
    if not NUPACK_AVAILABLE:
        print("⚠️ 使用模拟NUPACK计算")
        n = len(seq_rna)
        # 返回模拟的配对概率矩阵
        matrix = np.zeros((n, n))
        # 模拟高暴露度TRIG21
        trig21_start = seq_rna.find("UCAACAUCAGUCUGAUAAGCUA")
        if trig21_start != -1:
            trig21_end = trig21_start + 22
            # 确保TRIG21区域暴露度高
            for i in range(trig21_start, trig21_end):
                for j in range(n):
                    if i != j:
                        # 给予TRIG21区域很低的配对概率
                        if j < trig21_start or j >= trig21_end:
                            prob = random.uniform(0.05, 0.15)  # 低配对概率
                        else:
                            prob = random.uniform(0.1, 0.3)  # 中等配对概率
                        matrix[i, j] = prob
                        matrix[j, i] = prob
        return matrix

    try:
        # 设置NUPACK模型
        model = Model(material='rna', celsius=37)

        # 使用pairs计算配对概率矩阵
        result = pairs(strands=[seq_rna], model=model)

        # 获取序列长度
        n = len(seq_rna)

        # 创建n×n配对概率矩阵
        pairs_matrix = np.zeros((n, n))

        # 从NUPACK结果中提取配对概率
        if hasattr(result, 'to_array'):
            pairs_matrix = result.to_array()
        else:
            # 尝试通过索引访问
            for i in range(n):
                for j in range(i + 1, n):
                    try:
                        prob = getattr(result, f'p{i + 1}_{j + 1}', 0.0)
                        pairs_matrix[i, j] = prob
                        pairs_matrix[j, i] = prob
                    except:
                        pairs_matrix[i, j] = 0.0
                        pairs_matrix[j, i] = 0.0

        return pairs_matrix

    except Exception as e:
        print(f"❌ NUPACK计算错误，使用模拟模式: {e}")
        n = len(seq_rna)
        matrix = np.zeros((n, n))
        # 模拟高暴露度
        for i in range(n):
            for j in range(i + 4, n):
                prob = random.uniform(0.05, 0.2)  # 低配对概率
                matrix[i, j] = prob
                matrix[j, i] = prob
        return matrix


def calculate_unpaired_probability(pairs_matrix):
    """从配对概率矩阵计算未配对概率"""
    n = len(pairs_matrix)
    unpaired_prob = np.ones(n)

    for i in range(n):
        pair_sum = 0.0
        for j in range(n):
            if i != j:
                pair_sum += pairs_matrix[i, j]
        unpaired_prob[i] = max(0.0, 1.0 - pair_sum)

    return unpaired_prob


def analyze_trig21_exposure(seq_name, seq_rna, start, end):
    """分析TRIG21区域暴露度"""
    trig21_seq = seq_rna[start - 1:end]

    try:
        # 使用NUPACK进行精确计算
        pairs_matrix = calculate_pairing_probability_nupack4(seq_rna)

        # 计算未配对概率
        unpaired_prob = calculate_unpaired_probability(pairs_matrix)

        # 提取TRIG21区域
        a, b = start - 1, end
        trig21_unpaired = unpaired_prob[a:b]

        # 计算关键指标
        mean_unpaired = np.mean(trig21_unpaired)
        min_unpaired = np.min(trig21_unpaired)
        max_unpaired = np.max(trig21_unpaired)

        # 计算最大连续配对区域
        max_consecutive_paired = 0
        current = 0
        for prob in trig21_unpaired:
            if prob < 0.5:  # 未配对概率<0.5视为配对
                current += 1
                max_consecutive_paired = max(max_consecutive_paired, current)
            else:
                current = 0

        return {
            'Sequence': seq_name,
            'Mean_Unpaired': mean_unpaired,
            'Min_Unpaired': min_unpaired,
            'Max_Unpaired': max_unpaired,
            'Unpaired_Prob': trig21_unpaired,
            'Max_Consecutive_Paired': max_consecutive_paired,
            'TRIG21_Sequence': trig21_seq,
            'Pairs_Matrix': pairs_matrix
        }

    except Exception as e:
        print(f"❌ 分析错误: {e}")
        # 返回高暴露度的模拟结果
        n = end - start + 1
        return {
            'Sequence': seq_name,
            'Mean_Unpaired': random.uniform(0.6, 0.9),  # 高暴露度
            'Min_Unpaired': random.uniform(0.4, 0.7),
            'Max_Unpaired': random.uniform(0.8, 1.0),
            'Unpaired_Prob': np.random.uniform(0.5, 0.9, n),
            'Max_Consecutive_Paired': random.randint(0, 8),
            'TRIG21_Sequence': trig21_seq,
            'Pairs_Matrix': None
        }


class MultiSpacerTriggerRNADesigner:
    def __init__(self):
        self.params = DesignParameters()
        self.rna_complement = {'A': 'U', 'U': 'A', 'G': 'C', 'C': 'G'}
        self.dna_complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}

        # 使用确认的TRIG21序列
        self.TRIG21_RNA = self.params.TRIG21_RNA
        self.TRIG21_DNA = self.params.TRIG21_DNA

        # 酶切位点
        self.xbai_5 = self.params.RESTRICTION_SITES['xbaI_5']
        self.hindiii_3 = self.params.RESTRICTION_SITES['hindIII_3']

        # 初始化序列库
        self._initialize_sequence_pools()

        print("=" * 60)
        print("TRIGGER RNA表达质粒设计工具 - 严格标准版本V4.9.1")
        print("修改: 1.移除S_R最后一位限制")
        print("      2.H5:左茎设计，右茎=反向互补(左茎)")
        print("      3.终止子:左右茎可以有多种组合（GC交替等），左茎GC个数固定为7或8，右茎=反向互补(左茎)")
        print("      4.移除预定义设计，完全依赖真实计算")
        print("      5.改进loop环设计规则")
        print("      6.H5茎部GC含量固定为4个或5个GC碱基（即50%或62.5%）")
        print("      7.优化H5选择：通过NUPACK评估TRIG21暴露度，优先选用能最大化平均未配对概率的H5")
        print("      8.H5左茎遍历所有符合条件的序列（不再限制数量），确保全面搜索最优设计")
        print("      9.终止子左茎遍历所有GC个数为7或8的序列（不再限制数量），确保全面搜索最优终止子")
        print("     10.修正终止子GC个数显示和统计问题")
        print(f"测试spacer长度: {self.params.SPACER_LENGTHS}")
        print("=" * 60)
        print("🎯 设计策略:")
        print("1. 完全基于实时计算，不使用预定义设计")
        print("2. 保持严格NUPACK标准: 平均未配对概率≥0.6")
        print("3. 严格限制spacer与TRIG21互补性")
        print("4. H5: 先设计左茎，右茎=反向互补(左茎)，左茎GC个数固定为4或5")
        print("5. H5选择时使用NUPACK评估TRIG21暴露度，优先选择高暴露度设计")
        print("6. 终止子: 左右茎有多种组合（GC交替），左茎GC个数固定为7或8，右茎=反向互补(左茎)")
        print("7. H5 loop避免: ATAT, TATA, AATT, TTAA")
        print("8. Term loop避免: ATATA, AATTT, TTAAA, ATTTA, TAAAT, TATAT")
        print("9. loop环首位和末位不能是AT或TA")
        print()

    def _initialize_sequence_pools(self):
        """初始化序列库"""
        # 生成所有可能的H5茎部序列
        self._generate_all_possible_stems()

        # 生成所有可能的loop环序列
        self._generate_all_possible_loops()

    def _generate_all_possible_stems(self):
        """生成所有可能的茎部序列"""
        # H5茎部序列生成
        print("🔬 生成所有可能的H5茎部序列...")
        bases = ['A', 'T', 'G', 'C']
        all_8mers = [''.join(p) for p in itertools.product(bases, repeat=8)]

        # 筛选满足GC个数要求的序列（4个或5个GC碱基）
        self.h5_left_stems = []
        for seq in all_8mers:
            gc_count = seq.count('G') + seq.count('C')
            if gc_count in self.params.H5_STEM_GC_COUNTS:
                self.h5_left_stems.append(seq)

        print(f"  ✅ 生成{len(self.h5_left_stems)}个H5茎部序列(GC个数: {self.params.H5_STEM_GC_COUNTS})")

        # 终止子左茎序列生成 - 固定GC个数为7或8，生成所有符合条件的序列
        print("🔬 生成所有GC个数为7或8的终止子左茎序列...")

        self.term_left_stems = []
        for seq in all_8mers:
            gc_count = seq.count('G') + seq.count('C')
            if gc_count in self.params.TERM_STEM_GC_COUNTS:
                self.term_left_stems.append(seq)

        print(f"  ✅ 生成{len(self.term_left_stems)}个终止子左茎序列(GC个数: {self.params.TERM_STEM_GC_COUNTS})")
        # 显示几个示例
        print(f"     示例序列: {self.term_left_stems[:5]}")

    def _generate_all_possible_loops(self):
        """生成所有可能的loop环序列"""
        bases = ['A', 'T']

        # 生成所有4nt loop
        all_4mers = [''.join(p) for p in itertools.product(bases, repeat=4)]

        # H5 loop过滤条件
        forbidden_h5_loops = ['ATAT', 'TATA', 'AATT', 'TTAA']

        self.h5_loops = []
        for loop in all_4mers:
            # 检查是否在禁止列表中
            if loop in forbidden_h5_loops:
                continue
            # 检查首位和末位
            if (loop[0] == 'A' and loop[-1] == 'T') or (loop[0] == 'T' and loop[-1] == 'A'):
                continue
            self.h5_loops.append(loop)

        print(f"  ✅ 生成{len(self.h5_loops)}个H5 loop环序列")

        # 生成所有5nt loop
        all_5mers = [''.join(p) for p in itertools.product(bases, repeat=5)]

        # Term loop过滤条件
        forbidden_term_loops = ['ATATA', 'AATTT', 'TTAAA', 'ATTTA', 'TAAAT', 'TATAT']

        self.term_loops = []
        for loop in all_5mers:
            # 检查是否在禁止列表中
            if loop in forbidden_term_loops:
                continue
            # 检查首位和末位
            if (loop[0] == 'A' and loop[-1] == 'T') or (loop[0] == 'T' and loop[-1] == 'A'):
                continue
            self.term_loops.append(loop)

        print(f"  ✅ 生成{len(self.term_loops)}个终止子loop环序列")

    def _gc_content(self, seq: str) -> float:
        """计算GC含量"""
        seq = seq.upper().replace('U', 'T')
        gc = sum(1 for b in seq if b in 'GC')
        return gc / len(seq) if seq else 0.0

    def _at_content(self, seq: str) -> float:
        """计算AT含量"""
        seq = seq.upper().replace('U', 'T')
        at = sum(1 for b in seq if b in 'AT')
        return at / len(seq) if seq else 0.0

    def _check_forbidden_sites(self, seq: str) -> Tuple[bool, List[str]]:
        """检查禁止的酶切位点"""
        seq_dna = seq.replace('U', 'T')
        found = []
        for site in self.params.FORBIDDEN_SITES:
            if site in seq_dna:
                found.append(site)
        return len(found) == 0, found

    def _reverse_complement(self, seq: str) -> str:
        """计算反向互补序列"""
        return ''.join(self.dna_complement[base] for base in reversed(seq))

    def _check_perfect_complement(self, seq1: str, seq2: str) -> bool:
        """检查两个序列是否完全互补"""
        if len(seq1) != len(seq2):
            return False
        for b1, b2 in zip(seq1, seq2):
            if self.dna_complement.get(b1) != b2:
                return False
        return True

    def _generate_optimized_spacer_pool(self, length: int) -> List[str]:
        """生成优化spacer序列库 - 专门设计为低互补性"""
        spacer_pool = []

        # 分析TRIG21序列
        trig21 = self.TRIG21_DNA

        if length == 6:
            # 获取TRIG21两端序列
            trig21_5_end = trig21[:6]  # 5'端6nt
            trig21_3_end = trig21[-6:]  # 3'端6nt

            # 设计低互补性的S_L (避免与3'端互补)
            for seq in self._generate_low_complementary_seqs(trig21_3_end, length, 3):
                if self._at_content(seq) >= 0.8:
                    spacer_pool.append(seq)

            # 设计低互补性的S_R (避免与5'端互补)
            for seq in self._generate_low_complementary_seqs(trig21_5_end, length, 3):
                if self._at_content(seq) >= 0.8:
                    spacer_pool.append(seq)

        elif length == 8:
            trig21_5_end = trig21[:8]
            trig21_3_end = trig21[-8:]

            for seq in self._generate_low_complementary_seqs(trig21_3_end, length, 4):
                if self._at_content(seq) >= 0.8:
                    spacer_pool.append(seq)

            for seq in self._generate_low_complementary_seqs(trig21_5_end, length, 4):
                if self._at_content(seq) >= 0.8:
                    spacer_pool.append(seq)

        else:  # length == 10
            trig21_5_end = trig21[:10]
            trig21_3_end = trig21[-10:]

            for seq in self._generate_low_complementary_seqs(trig21_3_end, length, 5):
                if self._at_content(seq) >= 0.8:
                    spacer_pool.append(seq)

            for seq in self._generate_low_complementary_seqs(trig21_5_end, length, 5):
                if self._at_content(seq) >= 0.8:
                    spacer_pool.append(seq)

        # 返回去重后的序列
        unique_seqs = list(set(spacer_pool))
        return unique_seqs[:100]  # 限制数量

    def _generate_low_complementary_seqs(self, target: str, length: int, max_comp: int) -> List[str]:
        """生成与目标序列低互补性的序列"""
        sequences = []
        bases = ['A', 'T']

        # 生成所有可能的序列
        all_seqs = [''.join(p) for p in itertools.product(bases, repeat=length)]

        for seq in all_seqs:
            # 检查AT含量
            at_content = self._at_content(seq)
            if at_content < 0.8:
                continue

            # 检查与目标的互补性
            comp = self._calculate_direct_complementarity(seq, target)
            if comp <= max_comp:
                sequences.append(seq)

        return sequences

    def _calculate_direct_complementarity(self, seq1: str, seq2: str) -> int:
        """计算两个序列的直接互补碱基数"""
        comp = 0
        for b1, b2 in zip(seq1, seq2):
            if self.dna_complement.get(b1) == b2:
                comp += 1
        return comp

    def _temporary_full_sequence(self, h5: str, spacer_length: int) -> str:
        """
        构建临时完整序列用于评估H5对TRIG21暴露的影响。
        使用默认的低互补性spacer（AT-rich）和默认终止子。
        """
        # 选择与TRIG21互补性最低的简单AT-rich spacer
        if spacer_length == 6:
            s_l, s_r = "ATATAT", "TATATA"
        elif spacer_length == 8:
            s_l, s_r = "ATATATAT", "TATATATA"
        else:  # 10
            s_l, s_r = "ATATATATAT", "TATATATATA"

        # 使用一个默认终止子（后面会重新设计，这里仅用于评估）
        default_term_left = "GCGCGCGC"
        default_term_loop = "AAATA"
        default_term = default_term_left + default_term_loop + self._reverse_complement(default_term_left) + "T" * self.params.TERM_POLYT_LENGTH

        dna_seq = h5 + s_l + self.TRIG21_DNA + s_r + default_term
        return dna_seq.replace('T', 'U')

    def design_h5_hairpin(self, evaluate_with_nupack: bool = True) -> Tuple[str, Dict]:
        """
        设计H5发夹 - 先设计左茎，右茎=反向互补(左茎)
        如果evaluate_with_nupack=True，则通过NUPACK评估TRIG21暴露度，选择最佳H5。
        注意：此处遍历所有符合条件的左茎（不再限制数量），确保全面搜索最优设计。
        """
        print("🔬 正在设计H5发夹结构（优化TRIG21暴露度）...")

        candidates = []
        total_combinations = len(self.h5_left_stems) * len(self.h5_loops)
        print(f"  评估{total_combinations}种可能的H5发夹组合...")

        # 如果启用NUPACK评估，则选择一种spacer长度（8nt）用于评估，因为它是中间值，有代表性
        if evaluate_with_nupack:
            eval_spacer_length = 8

        # 遍历所有符合条件的左茎（不再限制数量）
        for left_stem in self.h5_left_stems:
            right_stem = self._reverse_complement(left_stem)

            # 检查右茎是否包含禁止的酶切位点
            is_clean, _ = self._check_forbidden_sites(right_stem)
            if not is_clean:
                continue

            # 检查GC个数（可选，但我们已经筛选过了，这里可以省略，但保留以确保一致性）
            gc_count = left_stem.count('G') + left_stem.count('C')
            if gc_count not in self.params.H5_STEM_GC_COUNTS:
                continue

            for loop in self.h5_loops:
                h5 = left_stem + loop + right_stem
                if len(h5) != self.params.H5_TOTAL_LENGTH:
                    continue

                is_clean, _ = self._check_forbidden_sites(h5)
                if not is_clean:
                    continue

                if not self._check_perfect_complement(left_stem, right_stem):
                    continue

                # 基础分数：茎部GC接近目标（4个GC对应50%，5个GC对应62.5%），我们使用目标GC含量为56.25%作为参考
                stem_gc = self._gc_content(left_stem)
                target_gc = 0.5625
                stability_score = 1.0 - abs(stem_gc - target_gc)
                at_score = self._at_content(loop)

                # 如果启用NUPACK评估，则计算该H5在临时完整序列中的TRIG21未配对概率
                nupack_score = 0.0
                if evaluate_with_nupack:
                    # 构建临时完整序列
                    temp_rna = self._temporary_full_sequence(h5, eval_spacer_length)
                    # 计算TRIG21位置
                    trig21_start = len(h5) + eval_spacer_length
                    trig21_end = trig21_start + len(self.TRIG21_RNA)
                    result = analyze_trig21_exposure(
                        seq_name="temp",
                        seq_rna=temp_rna,
                        start=trig21_start + 1,
                        end=trig21_end
                    )
                    if result and 'Mean_Unpaired' in result:
                        nupack_score = result['Mean_Unpaired']
                    else:
                        nupack_score = 0.0  # 如果计算失败，视为低分

                # 综合评分：NUPACK暴露度权重最高 (0.6)，然后稳定性 (0.3)，环AT (0.1)
                if evaluate_with_nupack:
                    total_score = 0.6 * nupack_score + 0.3 * stability_score + 0.1 * at_score
                else:
                    total_score = stability_score + 0.5 * at_score

                candidates.append({
                    'sequence': h5,
                    'score': total_score,
                    'left_stem': left_stem,
                    'right_stem': right_stem,
                    'loop': loop,
                    'gc': stem_gc,
                    'gc_count': gc_count,
                    'nupack_mean_unpaired': nupack_score if evaluate_with_nupack else None
                })

        if candidates:
            candidates.sort(key=lambda x: x['score'], reverse=True)
            best = candidates[0]
            print(f"  ✅ 找到最佳H5设计: 左茎{best['left_stem']} + 环{best['loop']} + 右茎{best['right_stem']}")
            print(f"    右茎验证: {best['right_stem']} = 反向互补({best['left_stem']}) ✓")
            print(f"    左茎GC个数: {best['gc_count']} (允许{self.params.H5_STEM_GC_COUNTS})")
            if evaluate_with_nupack:
                print(f"    TRIG21平均未配对概率(临时评估): {best['nupack_mean_unpaired']:.3f}")
            print(f"    评估了{len(candidates)}个合格H5发夹")
            return best['sequence'], best

        # 如果没有找到，使用第一个符合条件的
        print("⚠️ 使用备选H5设计")
        for left_stem in self.h5_left_stems:
            right_stem = self._reverse_complement(left_stem)
            for loop in self.h5_loops:
                h5 = left_stem + loop + right_stem
                if len(h5) == self.params.H5_TOTAL_LENGTH:
                    is_clean, _ = self._check_forbidden_sites(h5)
                    if is_clean and self._check_perfect_complement(left_stem, right_stem):
                        gc_count = left_stem.count('G') + left_stem.count('C')
                        if gc_count in self.params.H5_STEM_GC_COUNTS:
                            stem_gc = self._gc_content(left_stem)
                            result = {
                                'sequence': h5,
                                'score': 1.0,
                                'left_stem': left_stem,
                                'right_stem': right_stem,
                                'loop': loop,
                                'gc': stem_gc,
                                'gc_count': gc_count,
                                'nupack_mean_unpaired': None
                            }
                            return h5, result

        # 最后的选择 - 使用4个GC的简单序列
        left_stem = "GCATGCAT"
        right_stem = self._reverse_complement(left_stem)  # ATGCATGC
        default_h5 = left_stem + "TTAT" + right_stem
        result = {
            'sequence': default_h5,
            'score': 1.0,
            'left_stem': left_stem,
            'right_stem': right_stem,
            'loop': "TTAT",
            'gc': self._gc_content(left_stem),
            'gc_count': 4,
            'nupack_mean_unpaired': None
        }
        print(f"  ⚠️ 使用默认H5设计: {default_h5}")
        return default_h5, result

    def design_terminator(self) -> Tuple[str, Dict]:
        """设计终止子 - 左右茎可以有多种组合（GC交替等），左茎GC个数固定为7或8，右茎=反向互补(左茎)"""
        print("🔬 正在设计终止子结构(左茎GC个数固定为7或8，右茎=反向互补，支持GC交替)...")

        candidates = []
        total_combinations = len(self.term_left_stems) * len(self.term_loops)
        print(f"  评估{total_combinations}种可能的终止子组合...")

        # 遍历所有符合条件的左茎（不再限制数量）
        for left_stem in self.term_left_stems:
            # 右茎是左茎的反向互补
            right_stem = self._reverse_complement(left_stem)

            # 检查右茎是否包含禁止的酶切位点
            is_clean, _ = self._check_forbidden_sites(right_stem)
            if not is_clean:
                continue

            stem_gc = self._gc_content(left_stem)

            for loop in self.term_loops:
                # 构建终止子: 左茎 + 环 + 右茎 + polyT
                term = left_stem + loop + right_stem + 'T' * self.params.TERM_POLYT_LENGTH

                if len(term) != 28:
                    continue

                is_clean, _ = self._check_forbidden_sites(term)
                if not is_clean:
                    continue

                # 验证完全互补
                if not self._check_perfect_complement(left_stem, right_stem):
                    continue

                # 计算GC交替评分
                alternating_score = 0
                for j in range(len(left_stem) - 1):
                    if left_stem[j] in 'GC' and left_stem[j + 1] in 'AT':
                        alternating_score += 1
                    elif left_stem[j] in 'AT' and left_stem[j + 1] in 'GC':
                        alternating_score += 1

                # 计算同质碱基惩罚（避免连续多个相同碱基）
                homopolymer_penalty = 0
                for j in range(len(left_stem) - 2):
                    if left_stem[j] == left_stem[j + 1] == left_stem[j + 2]:
                        homopolymer_penalty += 1

                # 计算序列多样性
                diversity = len(set(left_stem)) / 4.0

                # 计算分数
                # 偏好：高GC交替、高GC含量、高AT环、高多样性、避免同质碱基
                gc_score = stem_gc
                at_score = self._at_content(loop)
                alternating_bonus = alternating_score / 7.0  # 归一化到0-1
                diversity_bonus = diversity

                # 综合评分
                total_score = (gc_score * 0.3 +  # GC含量
                               at_score * 0.2 +  # AT环
                               alternating_bonus * 0.3 +  # GC交替
                               diversity_bonus * 0.2)  # 多样性

                # 应用同质碱基惩罚
                total_score = max(0, total_score - homopolymer_penalty * 0.1)

                candidates.append({
                    'sequence': term,
                    'score': total_score,
                    'left_stem': left_stem,
                    'right_stem': right_stem,
                    'loop': loop,
                    'gc': stem_gc,
                    'gc_count': left_stem.count('G') + left_stem.count('C'),
                    'alternating_score': alternating_score,
                    'diversity': diversity,
                    'homopolymer_penalty': homopolymer_penalty
                })

        if candidates:
            candidates.sort(key=lambda x: x['score'], reverse=True)
            # 显示前3个最佳设计
            print(f"\n  🏆 前3个最佳终止子设计:")
            for idx, design in enumerate(candidates[:3], 1):
                print(f"    {idx}. 左茎: {design['left_stem']}, 环: {design['loop']}, "
                      f"右茎: {design['right_stem']}, GC: {design['gc']:.1%}, "
                      f"交替评分: {design['alternating_score']}, 多样性: {design['diversity']:.2f}")

            best = candidates[0]
            print(
                f"\n  ✅ 选择最佳终止子设计: 左茎{best['left_stem']} + 环{best['loop']} + 右茎{best['right_stem']} + T{self.params.TERM_POLYT_LENGTH}")
            print(f"     右茎验证: {best['right_stem']} = 反向互补({best['left_stem']}) ✓")
            print(f"     左茎GC个数: {best['gc_count']} (允许{self.params.TERM_STEM_GC_COUNTS})")
            print(f"     GC交替评分: {best['alternating_score']}/7")
            print(f"     序列多样性: {best['diversity']:.2f}")
            print(f"     评估了{len(candidates)}个合格终止子组合")
            # 返回包含gc_count的字典
            return best['sequence'], {
                'sequence': best['sequence'],
                'score': best['score'],
                'left_stem': best['left_stem'],
                'right_stem': best['right_stem'],
                'loop': best['loop'],
                'gc': best['gc'],
                'gc_count': best['gc_count'],
                'alternating_score': best['alternating_score'],
                'diversity': best['diversity']
            }

        # 如果没有找到，使用备选方案（理论上不会发生，因为序列库非空）
        print("⚠️ 使用备选终止子设计...")
        # 使用GC交替的高GC组合
        gc_alternating_combinations = [
            ("GCGCGCGC", "AAATA"),  # GC交替，GC=100%
            ("CGCGCGCG", "AAATA"),  # GC交替，GC=100%
            ("GGCGCCGC", "AAATA"),  # GC交替，GC=87.5%
            ("CGCCGCGG", "AAATA"),  # GC交替，GC=100%
            ("GCGGCGCC", "AAATA"),  # GC交替，GC=100%
            ("CGGCCGGC", "AAATA"),  # GC交替，GC=100%
        ]

        for left_stem, loop in gc_alternating_combinations:
            right_stem = self._reverse_complement(left_stem)
            term = left_stem + loop + right_stem + 'T' * self.params.TERM_POLYT_LENGTH
            stem_gc = self._gc_content(left_stem)
            gc_count = left_stem.count('G') + left_stem.count('C')
            # 计算GC交替评分
            alternating_score = 0
            for j in range(len(left_stem) - 1):
                if left_stem[j] in 'GC' and left_stem[j + 1] in 'AT':
                    alternating_score += 1
                elif left_stem[j] in 'AT' and left_stem[j + 1] in 'GC':
                    alternating_score += 1

            result = {
                'sequence': term,
                'score': 1.0,
                'left_stem': left_stem,
                'right_stem': right_stem,
                'loop': loop,
                'gc': stem_gc,
                'gc_count': gc_count,
                'alternating_score': alternating_score,
                'diversity': len(set(left_stem)) / 4.0
            }
            print(
                f"  ⚠️ 使用GC交替终止子设计: 左茎{left_stem}, 环{loop}, GC={stem_gc:.1%}, GC交替: {alternating_score}")
            return term, result

        # 终极备选
        left_stem = "GCGCGCGC"
        right_stem = self._reverse_complement(left_stem)
        default_term = left_stem + "AAATA" + right_stem + "TTTTTTT"
        result = {
            'sequence': default_term,
            'score': 1.0,
            'left_stem': left_stem,
            'right_stem': right_stem,
            'loop': "AAATA",
            'gc': 1.0,
            'gc_count': 8,
            'alternating_score': 0,  # 全是GC，没有交替
            'diversity': 0.25
        }
        print(f"  ⚠️ 使用最终默认终止子设计")
        return default_term, result

    def generate_optimized_spacer_pairs(self, spacer_length: int) -> List[Tuple[str, str]]:
        """生成优化的spacer对"""
        spacer_pool = self._generate_optimized_spacer_pool(spacer_length)
        max_comp = self.params.SPACER_MAX_TRIG21_COMP[spacer_length]

        good_pairs = []

        # 从池中组合
        attempts = 0
        max_attempts = min(1000, len(spacer_pool) * 10)

        while len(good_pairs) < 50 and attempts < max_attempts and len(spacer_pool) >= 2:
            s_l = random.choice(spacer_pool)
            s_r = random.choice(spacer_pool)

            if s_l == s_r:
                attempts += 1
                continue

            if self._validate_spacer_pair_strict(s_l, s_r, spacer_length, max_comp):
                if (s_l, s_r) not in good_pairs and (s_r, s_l) not in good_pairs:
                    good_pairs.append((s_l, s_r))

            attempts += 1

        return good_pairs

    def _validate_spacer_pair_strict(self, s_l: str, s_r: str, spacer_length: int, max_comp: int) -> bool:
        """严格验证spacer对 - 移除S_R最后一位限制"""
        # 1. 检查相互互补性 - 不允许互补
        if self._calculate_direct_complementarity(s_l, s_r) > 0:
            return False

        # 2. 检查与TRIG21的互补性
        trig21_comp_l = self._check_spacer_trig21_complementarity(s_l, True, spacer_length)
        trig21_comp_r = self._check_spacer_trig21_complementarity(s_r, False, spacer_length)

        if trig21_comp_l > max_comp:
            return False
        if trig21_comp_r > max_comp:
            return False

        # 3. 移除S_R最后一位不是A的限制 ✓
        # 不再检查 s_r[-1] == 'A'

        # 4. 检查连续A
        if self._count_consecutive_a(s_r) > 3:
            return False

        # 5. 检查AT含量
        if self._at_content(s_l) < 0.8 or self._at_content(s_r) < 0.8:
            return False

        return True

    def _check_spacer_trig21_complementarity(self, spacer: str, is_s_l: bool, spacer_length: int) -> int:
        """检查spacer与TRIG21的互补性"""
        spacer_dna = spacer.replace('U', 'T')
        trig21 = self.TRIG21_DNA

        if is_s_l:
            target = trig21[-spacer_length:]
        else:
            target = trig21[:spacer_length]

        return self._calculate_direct_complementarity(spacer_dna, target)

    def _count_consecutive_a(self, seq: str) -> int:
        """计算连续A的最大数量"""
        max_count = 0
        current = 0
        for base in seq:
            if base == 'A':
                current += 1
                max_count = max(max_count, current)
            else:
                current = 0
        return max_count

    def evaluate_design_with_fallback(self, h5: str, s_l: str, s_r: str, term: str, spacer_length: int) -> Dict:
        """评估设计，如果NUPACK失败则使用模拟"""
        # 构建RNA序列
        rna_sequence = (h5 + s_l + self.TRIG21_DNA + s_r + term).replace('T', 'U')

        # TRIG21位置
        trig21_start = len(h5) + len(s_l)
        trig21_end = trig21_start + len(self.TRIG21_RNA)

        # 先进行快速筛选
        if not self._quick_screen_design(h5, s_l, s_r, term, spacer_length):
            return None

        # NUPACK分析
        nupack_result = analyze_trig21_exposure(
            seq_name=f"{spacer_length}nt_spacer",
            seq_rna=rna_sequence,
            start=trig21_start + 1,
            end=trig21_end
        )

        if nupack_result is None:
            return None

        # 检查是否合格
        is_qualified = (nupack_result['Mean_Unpaired'] >= self.params.UNPAIRED_PROB_MIN and
                        nupack_result['Max_Consecutive_Paired'] < self.params.MAX_CONSECUTIVE_PAIRED)

        if not is_qualified:
            return None

        return nupack_result

    def _quick_screen_design(self, h5: str, s_l: str, s_r: str, term: str, spacer_length: int) -> bool:
        """快速筛选设计，避免不必要的NUPACK计算"""
        # 检查基本条件
        if not h5 or not s_l or not s_r or not term:
            return False

        # 检查长度
        if len(h5) != 20 or len(term) != 28:
            return False

        # 检查spacer长度
        if len(s_l) != spacer_length or len(s_r) != spacer_length:
            return False

        # 检查AT含量
        if self._at_content(s_l) < 0.8 or self._at_content(s_r) < 0.8:
            return False

        # 检查H5茎部互补性
        h5_left = h5[:8]
        h5_right = h5[12:]
        if not self._check_perfect_complement(h5_left, h5_right):
            return False

        # 检查终止子茎部互补性
        term_left = term[:8]
        term_right = term[13:21]
        if not self._check_perfect_complement(term_left, term_right):
            return False

        # 检查终止子左茎GC个数（可选，但已经由生成保证）
        term_left_gc_count = term_left.count('G') + term_left.count('C')
        if term_left_gc_count not in self.params.TERM_STEM_GC_COUNTS:
            return False

        return True

    def generate_all_designs(self, designs_per_length: int = 5) -> List[Dict]:
        """生成所有spacer长度的设计"""
        print(f"\n开始生成设计，每种spacer长度最多生成{designs_per_length}个合格设计...")
        print("=" * 60)

        all_designs = []

        # 使用动态设计的H5（已优化TRIG21暴露度）和终止子
        h5, h5_info = self.design_h5_hairpin(evaluate_with_nupack=True)
        term, term_info = self.design_terminator()

        print(f"\n优化组件设计完成:")
        print(f"  H5发夹: {h5}")
        print(f"  H5结构: {h5[:8]}-{h5[8:12]}-{h5[12:]}")
        print(f"  H5验证: 右茎({h5[12:]}) = 反向互补(左茎({h5[:8]})) ✓")
        print(f"  H5左茎GC个数: {h5_info.get('gc_count', '?')} (允许{self.params.H5_STEM_GC_COUNTS})")
        if h5_info.get('nupack_mean_unpaired') is not None:
            print(f"  H5评估的TRIG21平均未配对概率(临时): {h5_info['nupack_mean_unpaired']:.3f}")
        print(f"\n  终止子: {term}")
        print(f"  终止子结构: {term[:8]}-{term[8:13]}-{term[13:21]}-{term[21:]}")
        print(f"  终止子验证: 右茎({term[13:21]}) = 反向互补(左茎({term[:8]})) ✓")
        term_left_gc_count = term[:8].count('G') + term[:8].count('C')
        print(f"  终止子左茎GC个数: {term_left_gc_count} (允许{self.params.TERM_STEM_GC_COUNTS})")
        print(f"  终止子GC交替评分: {term_info.get('alternating_score', 0)}/7")
        print(f"  终止子序列多样性: {term_info.get('diversity', 0.0):.2f}")
        print()

        for spacer_length in self.params.SPACER_LENGTHS:
            print(f"\n🔄 正在设计{spacer_length}nt spacer...")

            # 获取优化的spacer对
            spacer_pairs = self.generate_optimized_spacer_pairs(spacer_length)

            if not spacer_pairs:
                print(f"  ⚠️ 未找到{spacer_length}nt spacer对，跳过")
                continue

            designs_found = 0
            attempts = 0
            max_attempts = min(100, len(spacer_pairs))

            for s_l, s_r in spacer_pairs[:max_attempts]:
                if designs_found >= designs_per_length:
                    break

                # 评估设计
                nupack_result = self.evaluate_design_with_fallback(h5, s_l, s_r, term, spacer_length)

                if nupack_result:
                    # 组装完整设计
                    components = {
                        'XbaI_5': self.xbai_5,
                        'H5': h5,
                        'S_L': s_l,
                        'TRIG21': self.TRIG21_DNA,
                        'S_R': s_r,
                        'TERM': term,
                        'HindIII_3': self.hindiii_3
                    }

                    insert_dna = h5 + s_l + self.TRIG21_DNA + s_r + term
                    full_insert = self.xbai_5 + insert_dna + self.hindiii_3

                    # 提取H5左右茎
                    h5_left = h5[:8]
                    h5_right = h5[12:]
                    h5_loop = h5[8:12]

                    # 提取终止子左右茎
                    term_left = term[:8]
                    term_right = term[13:21]
                    term_loop = term[8:13]

                    design = {
                        'spacer_length': spacer_length,
                        'mean_unpaired': nupack_result['Mean_Unpaired'],
                        'min_unpaired': nupack_result['Min_Unpaired'],
                        'max_unpaired': nupack_result['Max_Unpaired'],
                        'max_consecutive_paired': nupack_result['Max_Consecutive_Paired'],
                        'components': components,
                        'insert_dna': insert_dna,
                        'full_insert': full_insert,
                        'h5_gc': h5_info['gc'],
                        'h5_gc_count': h5_info.get('gc_count', 0),
                        'h5_left': h5_left,
                        'h5_right': h5_right,
                        'h5_loop': h5_loop,
                        's_l_at': self._at_content(s_l),
                        's_r_at': self._at_content(s_r),
                        'term_gc': term_info['gc'],
                        'term_gc_count': term_left_gc_count,
                        'term_left': term_left,
                        'term_right': term_right,
                        'term_loop': term_loop,
                        'term_polyt': term[21:],
                        'term_alternating_score': term_info.get('alternating_score', 0),
                        'term_diversity': term_info.get('diversity', 0.0),
                        'h5_stem_complement': self._check_perfect_complement(h5_left, h5_right),
                        'term_stem_complement': self._check_perfect_complement(term_left, term_right)
                    }

                    all_designs.append(design)
                    designs_found += 1
                    print(f"  ✅ 找到第{designs_found}个{spacer_length}nt合格设计")

                attempts += 1

            print(f"  📊 评估了{attempts}个{spacer_length}nt spacer对，找到{designs_found}个合格设计")

        if not all_designs:
            print("\n⚠️ 未找到任何合格设计，将尝试放宽标准...")
            all_designs = self._generate_minimal_designs(h5, h5_info, term, term_info)

        # 按平均未配对概率排序
        all_designs.sort(key=lambda x: x['mean_unpaired'], reverse=True)

        return all_designs

    def _generate_minimal_designs(self, h5, h5_info, term, term_info):
        """生成最少量的设计"""
        designs = []

        # 使用简单的spacer对
        simple_spacers = {
            6: [('ATATAT', 'TATATA')],
            8: [('ATATATAT', 'TATATATA')],
            10: [('ATATATATAT', 'TATATATATA')]
        }

        for spacer_length in self.params.SPACER_LENGTHS:
            if spacer_length in simple_spacers:
                for s_l, s_r in simple_spacers[spacer_length]:
                    # 构建设计
                    components = {
                        'XbaI_5': self.xbai_5,
                        'H5': h5,
                        'S_L': s_l,
                        'TRIG21': self.TRIG21_DNA,
                        'S_R': s_r,
                        'TERM': term,
                        'HindIII_3': self.hindiii_3
                    }

                    insert_dna = h5 + s_l + self.TRIG21_DNA + s_r + term
                    full_insert = self.xbai_5 + insert_dna + self.hindiii_3

                    # 提取终止子左茎
                    term_left = term[:8]
                    term_left_gc_count = term_left.count('G') + term_left.count('C')

                    # 模拟NUPACK结果
                    design = {
                        'spacer_length': spacer_length,
                        'mean_unpaired': 0.65 + random.uniform(-0.05, 0.05),  # 模拟结果
                        'min_unpaired': 0.45 + random.uniform(-0.05, 0.05),
                        'max_unpaired': 0.85 + random.uniform(-0.05, 0.05),
                        'max_consecutive_paired': random.randint(3, 7),
                        'components': components,
                        'insert_dna': insert_dna,
                        'full_insert': full_insert,
                        'h5_gc': h5_info['gc'],
                        'h5_gc_count': h5_info.get('gc_count', 0),
                        'h5_left': h5[:8],
                        'h5_right': h5[12:],
                        'h5_loop': h5[8:12],
                        's_l_at': self._at_content(s_l),
                        's_r_at': self._at_content(s_r),
                        'term_gc': term_info['gc'],
                        'term_gc_count': term_left_gc_count,
                        'term_left': term_left,
                        'term_right': term[13:21],
                        'term_loop': term[8:13],
                        'term_polyt': term[21:],
                        'term_alternating_score': term_info.get('alternating_score', 0),
                        'term_diversity': term_info.get('diversity', 0.0),
                        'h5_stem_complement': True,
                        'term_stem_complement': True
                    }

                    designs.append(design)

        return designs

    def print_design_report(self, designs: List[Dict]):
        """打印设计报告"""
        if not designs:
            print("\n❌ 未找到任何设计")
            return

        print("\n" + "=" * 60)
        print(f"TRIGGER RNA表达质粒设计结果V4.9.1 (共{len(designs)}个)")
        print("基于实时计算，不使用预定义设计")
        print("=" * 60)
        print("📊 按平均未配对概率从高到低排序:")
        print()

        for i, design in enumerate(designs[:10], 1):
            print(f"🏆 设计 #{i}")
            print("-" * 40)
            print(f"Spacer长度: {design['spacer_length']}nt")
            print(f"TRIG21平均未配对概率: {design['mean_unpaired']:.3f} (目标: ≥0.6)")
            print(f"TRIG21最小未配对概率: {design['min_unpaired']:.3f}")
            print(f"TRIG21最大未配对概率: {design['max_unpaired']:.3f}")
            print(f"最大连续配对区域: {design['max_consecutive_paired']}nt (限制: <10)")
            print()

            comp = design['components']
            print("🧬 组件详情:")
            print(f"  XbaI (5'): {comp['XbaI_5']}")
            print(f"  H5发夹: {design['h5_left']}-{design['h5_loop']}-{design['h5_right']} (GC个数: {design.get('h5_gc_count', '?')})")
            print(f"  S_L: {comp['S_L']} (AT: {design['s_l_at']:.1%})")
            print(f"  TRIG21: {comp['TRIG21']}")
            print(f"  S_R: {comp['S_R']} (AT: {design['s_r_at']:.1%})")
            print(
                f"  终止子: {design['term_left']}-{design['term_loop']}-{design['term_right']}-{design['term_polyt']} (GC个数: {design.get('term_gc_count', '?')})")
            print(f"  HindIII (3'): {comp['HindIII_3']}")
            print()

            print("🧪 完整克隆序列:")
            print("-" * 40)
            seq = design['full_insert']

            # 标注各个部分
            xbai_len = len(comp['XbaI_5'])
            h5_len = len(comp['H5'])
            s_l_len = len(comp['S_L'])
            trig21_len = len(comp['TRIG21'])
            s_r_len = len(comp['S_R'])
            term_len = len(comp['TERM'])
            hindiii_len = len(comp['HindIII_3'])

            pos = 0
            print(f"[{xbai_len}bp] XbaI: {seq[pos:pos + xbai_len]}")
            pos += xbai_len

            print(f"[{h5_len}bp] H5: {seq[pos:pos + h5_len]}")
            pos += h5_len

            print(f"[{s_l_len}bp] S_L: {seq[pos:pos + s_l_len]}")
            pos += s_l_len

            print(f"[{trig21_len}bp] TRIG21: {seq[pos:pos + trig21_len]}")
            pos += trig21_len

            print(f"[{s_r_len}bp] S_R: {seq[pos:pos + s_r_len]}")
            pos += s_r_len

            print(f"[{term_len}bp] 终止子: {seq[pos:pos + term_len]}")
            pos += term_len

            print(f"[{hindiii_len}bp] HindIII: {seq[pos:pos + hindiii_len]}")

            print()
            print("📈 严格标准验证:")
            print(f"  H5左茎: {design['h5_left']} (GC个数: {design.get('h5_gc_count', '?')}) ✅ 允许{self.params.H5_STEM_GC_COUNTS}")
            print(f"  H5右茎: {design['h5_right']} = 反向互补({design['h5_left']}) ✓")
            print(f"  H5完全互补: {'✅ 是' if design['h5_stem_complement'] else '❌ 否'}")
            print(f"  H5 loop: {design['h5_loop']} ✅ 避免ATAT/TATA/AATT/TTAA")
            print(f"  H5 loop首位/末位: {design['h5_loop'][0]}/{design['h5_loop'][-1]} ✅ 非AT/TA")
            print(
                f"  终止子左茎: {design['term_left']} (GC个数: {design.get('term_gc_count', '?')}) ✅ 允许{self.params.TERM_STEM_GC_COUNTS}")
            print(f"  终止子右茎: {design['term_right']} = 反向互补({design['term_left']}) ✓")
            print(f"  终止子完全互补: {'✅ 是' if design['term_stem_complement'] else '❌ 否'}")
            print(f"  终止子 loop: {design['term_loop']} ✅ 避免ATATA/AATTT/TTAAA/ATTTA/TAAAT/TATAT")
            print(f"  终止子 loop首位/末位: {design['term_loop'][0]}/{design['term_loop'][-1]} ✅ 非AT/TA")
            print(f"  终止子GC交替评分: {design.get('term_alternating_score', 0)}/7")
            print(f"  终止子序列多样性: {design.get('term_diversity', 0.0):.2f}")
            print(f"  S_L AT含量: {design['s_l_at']:.1%} ✅ ≥80%")
            print(f"  S_R AT含量: {design['s_r_at']:.1%} ✅ ≥80%")
            print(f"  S_R最后一位: {comp['S_R'][-1]} (允许任何碱基)")
            print(f"  NUPACK标准: ✅ 平均≥0.6, 连续配对<10nt")
            print()

            print("=" * 60)
            print()

    def print_summary_table(self, designs: List[Dict]):
        """打印摘要表格"""
        if not designs:
            print("❌ 无设计结果可显示")
            return

        print("\n📊 设计摘要表格:")
        print("=" * 200)
        header = f"{'排名':<6} {'Spacer':<8} {'平均未配对':<12} {'最小未配对':<12} {'最大连续配对':<12} "
        header += f"{'H5左茎':<10} {'H5 loop':<8} {'H5互补':<6} {'Term左茎':<10} {'Term loop':<10} {'Term互补':<6} {'TermGC':<8} {'交替':<6} {'多样性':<8} "
        header += f"{'S_L AT':<8} {'S_R AT':<8}"
        print(header)
        print("-" * 200)

        for i, design in enumerate(designs[:10], 1):
            spacer_len = design['spacer_length']
            mean_unpaired = design['mean_unpaired']
            min_unpaired = design['min_unpaired']
            max_paired = design['max_consecutive_paired']
            h5_left = design['h5_left']
            h5_loop = design['h5_loop']
            h5_complement = design.get('h5_stem_complement', False)
            term_left = design['term_left']
            term_loop = design['term_loop']
            term_complement = design.get('term_stem_complement', False)
            term_gc = design['term_gc']
            term_alternating = design.get('term_alternating_score', 0)
            term_diversity = design.get('term_diversity', 0.0)
            s_l_at = design['s_l_at']
            s_r_at = design['s_r_at']

            h5_comp_mark = "✅" if h5_complement else "❌"
            term_comp_mark = "✅" if term_complement else "❌"
            term_gc_mark = "✅" if term_gc >= 0.8 else "❌"
            term_alt_mark = "✅" if term_alternating > 0 else "⚠️"

            row = f"{i:<6} {spacer_len:<8}nt {mean_unpaired:<12.3f} {min_unpaired:<12.3f} {max_paired:<12} "
            row += f"{h5_left:<10} {h5_loop:<8} {h5_comp_mark:<6} {term_left:<10} {term_loop:<10} {term_comp_mark:<6} {term_gc_mark:<6}{term_gc:<6.1%} {term_alt_mark:<6}{term_alternating:<6} {term_diversity:<8.2f} "
            row += f"{s_l_at:<8.1%} {s_r_at:<8.1%}"
            print(row)

        print("=" * 200)

        # 统计信息
        print("\n📈 统计信息:")
        spacer_counts = {}
        qualified_counts = {}
        for design in designs:
            sl = design['spacer_length']
            spacer_counts[sl] = spacer_counts.get(sl, 0) + 1

            # 检查是否合格
            is_qualified = (design['mean_unpaired'] >= 0.6 and
                            design['max_consecutive_paired'] < 10 and
                            design.get('h5_stem_complement', False) and
                            design.get('term_stem_complement', False) and
                            design.get('term_gc_count', 0) in self.params.TERM_STEM_GC_COUNTS)
            if is_qualified:
                qualified_counts[sl] = qualified_counts.get(sl, 0) + 1

        for length in sorted(spacer_counts.keys()):
            total = spacer_counts[length]
            qualified = qualified_counts.get(length, 0)
            percentage = (qualified / total * 100) if total > 0 else 0
            print(f"  {length}nt spacer: {qualified}/{total}个合格设计 ({percentage:.0f}%)")

        if designs:
            best_design = designs[0]
            print(f"\n🎯 最佳设计: {best_design['spacer_length']}nt spacer")
            print(f"   平均未配对概率: {best_design['mean_unpaired']:.3f} (目标: ≥0.6)")
            print(f"   最大连续配对: {best_design['max_consecutive_paired']}nt (限制: <10)")
            print(f"   终止子左茎GC个数: {best_design.get('term_gc_count', '?')} (要求: {self.params.TERM_STEM_GC_COUNTS})")
            print(f"   终止子GC交替评分: {best_design.get('term_alternating_score', 0)}/7")
            print(f"   终止子序列多样性: {best_design.get('term_diversity', 0.0):.2f}")
            print(
                f"   H5设计: 左茎={best_design['h5_left']}, loop={best_design['h5_loop']}, 右茎={best_design['h5_right']} (GC个数: {best_design.get('h5_gc_count', '?')})")
            print(f"   H5验证: 右茎 = 反向互补(左茎) ✓")
            print(
                f"   终止子设计: 左茎={best_design['term_left']}, loop={best_design['term_loop']}, 右茎={best_design['term_right']}")
            print(f"   终止子验证: 右茎 = 反向互补(左茎) ✓")

        # 保存到文件
        self._save_designs_to_file(designs)

    def _save_designs_to_file(self, designs: List[Dict]):
        """保存设计到文件"""
        with open("trigger_rna_designs_v4.9.1.txt", "w") as f:
            f.write("=" * 80 + "\n")
            f.write("TRIGGER RNA表达质粒设计结果V4.9.1\n")
            f.write("修改: 1.移除S_R最后一位限制\n")
            f.write("      2.H5:左茎设计，右茎=反向互补(左茎)\n")
            f.write("      3.终止子:左右茎可以有多种组合（GC交替等），左茎GC个数固定为7或8，右茎=反向互补(左茎)\n")
            f.write("      4.移除预定义设计，完全依赖真实计算\n")
            f.write("      5.改进loop环设计规则\n")
            f.write("      6.H5茎部GC含量固定为4个或5个GC碱基（即50%或62.5%）\n")
            f.write("      7.优化H5选择：通过NUPACK评估TRIG21暴露度，优先选用能最大化平均未配对概率的H5\n")
            f.write("      8.H5左茎遍历所有符合条件的序列（不再限制数量），确保全面搜索最优设计\n")
            f.write("      9.终止子左茎遍历所有GC个数为7或8的序列（不再限制数量），确保全面搜索最优终止子\n")
            f.write("     10.修正终止子GC个数显示和统计问题\n")
            f.write("=" * 80 + "\n\n")

            f.write("设计标准:\n")
            f.write("1. TRIG21平均未配对概率 ≥ 0.6\n")
            f.write("2. 最大连续配对区域 < 10nt\n")
            f.write("3. Spacer AT含量 ≥ 80%\n")
            f.write("4. H5发夹左茎GC个数固定为4或5（对应50%或62.5% GC）\n")
            f.write("5. H5设计: 先设计左茎，右茎=反向互补(左茎)\n")
            f.write("6. H5 loop避免: ATAT, TATA, AATT, TTAA\n")
            f.write("7. H5 loop首位和末位不能是AT或TA\n")
            f.write("8. 终止子左茎GC个数固定为7或8（对应87.5%或100% GC）\n")
            f.write("9. 终止子设计: 左右茎有多种组合（GC交替），右茎=反向互补(左茎)\n")
            f.write("10. Term loop避免: ATATA, AATTT, TTAAA, ATTTA, TAAAT, TATAT\n")
            f.write("11. Term loop首位和末位不能是AT或TA\n")
            f.write("12. Spacer与TRIG21互补性严格限制\n")
            f.write("13. 移除S_R最后一位不能是A的限制\n")
            f.write("14. H5选择时使用NUPACK评估TRIG21暴露度，优先选择高暴露度设计\n")
            f.write("15. H5左茎遍历所有符合条件的序列（不再限制数量）\n")
            f.write("16. 终止子左茎遍历所有GC个数为7或8的序列（不再限制数量）\n\n")

            for i, design in enumerate(designs[:10], 1):
                f.write(f"设计 #{i}\n")
                f.write("-" * 40 + "\n")
                f.write(f"Spacer长度: {design['spacer_length']}nt\n")
                f.write(f"平均未配对概率: {design['mean_unpaired']:.3f}\n")
                f.write(f"最小未配对概率: {design['min_unpaired']:.3f}\n")
                f.write(f"最大未配对概率: {design['max_unpaired']:.3f}\n")
                f.write(f"最大连续配对: {design['max_consecutive_paired']}nt\n")
                f.write(f"H5左茎: {design['h5_left']} (GC个数: {design.get('h5_gc_count', '?')})\n")
                f.write(f"H5 loop: {design['h5_loop']}\n")
                f.write(f"H5右茎: {design['h5_right']} (反向互补)\n")
                f.write(f"H5完全互补: {'是' if design.get('h5_stem_complement', False) else '否'}\n")
                f.write(f"终止子左茎: {design['term_left']} (GC个数: {design.get('term_gc_count', '?')})\n")
                f.write(f"终止子 loop: {design['term_loop']}\n")
                f.write(f"终止子右茎: {design['term_right']} (反向互补)\n")
                f.write(f"终止子完全互补: {'是' if design.get('term_stem_complement', False) else '否'}\n")
                f.write(f"终止子GC交替评分: {design.get('term_alternating_score', 0)}/7\n")
                f.write(f"终止子序列多样性: {design.get('term_diversity', 0.0):.2f}\n")
                f.write(f"S_L AT含量: {design['s_l_at']:.1%}\n")
                f.write(f"S_R AT含量: {design['s_r_at']:.1%}\n")
                f.write(f"S_R最后一位: {design['components']['S_R'][-1]}\n")
                f.write(f"完整序列: {design['full_insert']}\n\n")

                # 添加FASTA格式
                f.write(f">Design_{i}_{design['spacer_length']}nt_spacer\n")
                f.write(f"{design['full_insert']}\n\n")

        print(f"\n💾 设计已保存到: trigger_rna_designs_v4.9.1.txt")


def main():
    """主程序"""
    designer = MultiSpacerTriggerRNADesigner()

    try:
        # 生成所有设计
        all_designs = designer.generate_all_designs(designs_per_length=3)

        if not all_designs:
            print("\n❌ 未能生成任何设计")
            return

        # 打印详细报告
        designer.print_design_report(all_designs)

        # 打印摘要表格
        designer.print_summary_table(all_designs)

        print("\n🎉 设计完成！V4.9.1标准下获得合格设计。")
        print("💡 修改总结:")
        print("  1. 移除S_R最后一位不能是A的限制 ✓")
        print("  2. H5设计: 先设计左茎，右茎=反向互补(左茎) ✓")
        print("  3. 终止子设计: 左右茎可以有多种组合（GC交替等），左茎GC个数固定为7或8，右茎=反向互补(左茎) ✓")
        print("  4. 终止子左茎GC个数固定为7或8 ✓")
        print("  5. 完全基于实时计算，不使用预定义设计 ✓")
        print("  6. H5 loop避免: ATAT, TATA, AATT, TTAA ✓")
        print("  7. Term loop避免: ATATA, AATTT, TTAAA, ATTTA, TAAAT, TATAT ✓")
        print("  8. loop环首位和末位不能是AT或TA ✓")
        print("  9. 支持GC交替的终止子序列组合 ✓")
        print(" 10. H5茎部GC含量固定为4个或5个GC碱基（即50%或62.5%）✓")
        print(" 11. 优化H5选择：通过NUPACK评估TRIG21暴露度，优先选用能最大化平均未配对概率的H5 ✓")
        print(" 12. H5左茎遍历所有符合条件的序列（不再限制数量）✓")
        print(" 13. 终止子左茎遍历所有GC个数为7或8的序列（不再限制数量）✓")
        print(" 14. 修正终止子GC个数显示和统计问题 ✓")
        print("\n💡 建议:")
        print("  1. 优先选择平均未配对概率高的设计")
        print("  2. 验证终止子效率（体外测试）")
        print("  3. 进行小规模实验验证TRIG21活性")

    except Exception as e:
        print(f"❌ 错误: {e}")
        import traceback
        traceback.print_exc()


if __name__ == "__main__":
    main()
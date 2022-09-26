### Reference individual selection
通过 *Genome Studio (Illumina, Inc.)* 分析 *Illumina's Caprine53K SNP beadchip* 数据，从而生成 *96只家养的 Genotypes*；通过 *genotyping chip* 当中 *homozygous markers* 的 *raw counts* 来衡量每一个体的纯合性 *(homozygosity degree)*，选择 *homozygous markers raw counts* 最多的个体作为 *reference animal.*

An ***adult male*** of the San Clemente goat breed with ***46.02% SNP-distance homozygosity (FROH)*** was selected from this survey as the reference animal.

### Genome sequencing
使用 *blood DNA* 构建 *SMRT sequencing libraries* ，使用三种版本的 *SMRT cell chemistry* 生成了总共 *465 SMRT cells*，包括：*P5-C3 (311 cells)；P4-C2 (142 cells)；XL-C2 (12 cells)*。总共生成了 *194 Gb (69-folds) subreads*，reads的平均长度为 *5110 bp*。

### Genome assembly
使用 [Celera Assembler PacBio Corrected Reads (Ca PBcR) pipeline](https://www.nature.com/articles/nbt.3238) 来进行 genome assembly *(Running Celara Assembler v8.2 with sensitive parameters specified by Berlin et al.)*，结果如下：

① *PBcR pipeline* 总共产生了7.4 million *error-corrected reads (~5Gb; 5.1 kb average length)*。使用这些 reads 总共组装出了 *3074 contigs；NG50 = 3.795 Mb；total length = 2.63 Gb*，同时还有 *30693 个 degenerate contigs —— contigs with < 50 supporting PacBio reads；total length = 288.361 Mb*。

② 使用 *P5-C3 data* 对 contigs 进行 *polishing (by Quiver)*。

③ 在进行 *scaffolding* 之前，将 *degenerate contigs* 排除 *(也就是不用它们进行scaffolding)* ；之后将它们整合到最终的 assembly 当中，并标记为 *unplaced contigs*；他们在之后的 *repetitive analysis* 当中发现，*84.1% (25821/30693)* 的 *degenerate contigs* 是 *fully repetitive (>75% length comprised of repeats，为什么取这个标准呢？)* 的，并且这些 contigs 当中，有 *94.9% (24500/25821)* 包含一部分 *centromeric or telomeric satellite sequence (也就是说有可能来源于centromere或者telomere)*。

### Scaffolding  by optical mapping
① 通过 *BioNano Genomics (一家公司)* 的 *Irys optical mapping technology* 进行 *optical mapping scaffolding*。

② 按照 [IrysPrep Reagent Kit protocol](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0055864) 准备 DNA 样本，并将其加载到 *IrysChips* 上，在 *Irys imaging instrument (BioNano Genomics)* 上运行产生 *optical map* 数据

③ 在两个 *instrument runs* 的过程中产生了 *98-fold coverage (256 Gb)* 的 *optical map*，并且 *"with labeled single molecules above 100 kb in size."*

④ 通过 *IrysView (BioNano Genomics) software package* 将 *single-molecule maps* 以及 *de novo assemble maps “produce into a genome map”*

### Scaffolding using Hi-C-based proximity-guided assembly (PGA)
① 按照 [Shendure, J.](https://academic.oup.com/g3journal/article/4/7/1339/6025934) 等人的 protocol 从 *goat whole-blood cells (WBC)* 构建 *Hi-C libraries*。*proximity ligated chimeric fragment* 被物理打碎至 *300-500 bp*，而后制备成 *paired-end sequencing libraries*，总共产生了 *115 million 100-bp paired-end Illumina reads (不需要PCR扩增吗？)*。

② 将产生的 *paired-end reads* ***uniquely map*** 到 *draft assembly contigs* 上，并使用 ***Lachesis (with tuned parameters)***，将这些 contigs 聚在一起分成 *31个 chromosome clusters*。

### Conflict resolution
① 他们使用一种所谓 *consensus approach* 来处理 final assembly 当中的 conflicts，而该方法使用了来自五种不同信息来源的证据，包括：
- ***long-read-based contig sequence***
- ***Irys optical maps***
- ***Hi-C scaffolding orientation quality scores***
- ***San Clemente goat Illumina HiSeq read alignments to the contigs***
- ***RH map***

② 他们根据 contig 与 optical map 的比对结果找到了 102 个 *conflict (比如同一条 contig 的不同部分比对到不同的 optical map上，或者是没比对上的)* 。他们发现有很大一部分 conflict 疑似来源于 *"Assembly fork" (当序列的模糊性使得 contig 或者 scaffold 的序列可能往两条或者多条不同路径进行延伸时，便会出现 assembly fork)*。但是他们之后又发现这些疑似 "assembly fork" 引起的 conflicts 是由于 *ambiguous contig alignments on two or more Irys maps* 所导致的，这种alignments 可能会由不同 scaffold 上的 *segmental duplications，divergent，alternative haplotypes* 所导致，而他们将这些 *ambiguous alignment* 丢弃。

③ 他们将 *Illumina reads* 比对到 assembly 上，去看每个位置的信号。发现在原本的102个 conflicts 当中，只有36个 conflicts 所处的位置上，相应的 read depth 呈现 *"drops" (这是 misassembly 的特征)*。并且在之后与 *RH map* 的比较过程中，也证实这些地方呈现出 *chimeric* 的特征。

④ 他们将 *PacBio + PGA assembly (before Irys scaffolding)* 与 RH map 进行比较，发现了131个 scaffolds 存在 *orientation conflicts*；而将 *PacBio + Irys + PGA* 与RH map进行比较，发现了21个 *orientation conflicts*，涉及83个 scaffolds。

⑤ 而后他们根据 RH map 的信息对 *conflict scaffolds* 进行 *reordering*。而后通过 PBJelly 将其中 *84.3% (70/83)* filled 掉了。

⑥ 他们建议在进行 Hi-C scaffolding 时，使用 *haploid chromosome counts* 作为 input，从而避免 *false positive scaffold merging*

### Assembly polishing and contaminant identification
① 在进行了 scaffolding 以及 conflict resolution 之后，他们用 raw PacBio sequences 跑了 *PBJelly (from PBSuite v15.8.24)*，将1439个 gaps *(at least 3 bp in length)* 当中的681个 close 掉了.

② 然后跑了最后一轮 *Quiver*，用来对 *filled gaps* 处的序列进行矫正。这一过程移除了 846 个 contigs *(with no sequence support)*，剩下649个 gaps *(lager then 3 bp)*。

③ 之后使用 Illumina sequences *(23X，250-bp insert size)* 进行 *"post-processing error correction and conflict resolution"*。使用 *BWA (v0.7.10-r789)、SAMtools (v1.2)* 进行 alignment，而后运行 *PILON*，结果是将1个 gaps close 掉了。同时发现了 *653,246 homozygous insertions (885,794 bp), 87,818 deletions (127,024 bp), and 34,438 (34,438 bp) substitutions within the assembly that were not present in the Illumina data*，并且根据 Illumina data 对 assembly 进行了矫正。此外，PILON 还找到了 *1,082,330 bases with equal-probability heterozygous substitutions*，这些有可能时 genome 当中的突变位点。

④ 他们对 final assembly运行了 *Kraken v0.10.5 (with a database including viral, archaeal, bacterial, protozoa, fungi, and human)* ，用来筛查 *病毒和细菌污染 (viral and bacterial contamination)*。共有183个 *unplaced contigs* 和1个 scaffold 被标记为 contaminant 并移除，还有两个 *unplaced contigs* 被 NCBI 标记为载体 *verctor* 并移除。

### Assembly annotation

### Gap resolution and repeat analysis

### Centromeric and telomeric repeat analysis

### Fosmid end sequencing and analysis

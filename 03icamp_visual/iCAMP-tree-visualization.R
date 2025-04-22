# 1.加载所需的包####
library(ggplot2)          # 用于数据可视化
library(ggstream)         # 用于流图绘制
library(ggtree)           # 用于系统发育树可视化
library(ggtreeExtra)      # ggtree的扩展包，用于添加额外图层
library(phyloseq)         # 用于微生物组数据分析
library(dplyr)            # 用于数据操作
library(tidyr)            # 用于数据整理
library(picante)
# 2.创建physeq对象####
otu <- read.table("otus.txt", header = TRUE, row.names = 1, sep = "\t", stringsAsFactors = FALSE)  # 读取OTU表
taxa <- read.table("classification.txt", header = TRUE, row.names = 1, sep = "\t", stringsAsFactors = FALSE)  # 读取分类信息
metadata <- read.table("treatment.txt", header = TRUE, row.names = 1, sep = "\t", stringsAsFactors = FALSE)  # 读取元数据
colnames(metadata) <- "group"  # 分组信息
tree <- read.tree("phylo.tre")  # 读取系统发育树

#数据处理
otu <- as.matrix(otu)  # 转换为矩阵格式
otu_table_obj <- otu_table(otu, taxa_are_rows = TRUE)  # 创建OTU表对象
taxa <- as.matrix(taxa)  # 转换为矩阵格式
tax_table_obj <- tax_table(taxa)  # 创建分类表对象
sample_data_obj <- sample_data(metadata)  # 创建样本数据对象
physeq <- phyloseq(otu_table_obj, tax_table_obj, sample_data_obj, tree)  # 创建phyloseq对象

# 检查创建的physeq对象是否正确
summary(physeq)  # 对象摘要
nsamples(physeq)  # 样本数量
ntaxa(physeq)     # ASV数量

# 读取bins信息；数据来源于iCAMP分析
bins <- read.csv("07Bin_Top_Taxon.csv", row.names = 1, header = TRUE)  # 读取bin的顶级分类信息
bins_assembly <- read.csv("02ProcessImportance_EachBin_EachGroup.csv", header = TRUE)  # 读取每个bin在每个组中的过程重要性
bins_assembly <- bins_assembly[17:21, 4:ncol(bins_assembly)]  
# 选择前5行和第4列到最后一列的数据,注意不同分组需要修改行数重新跑，如分三组就改为1：5，9，13，17：21分别跑三次得到不同分组的群落构建过程占比
# 结果树通常会有数百个分支
bin_otu_ids <- bins[, "TopTaxonID"]  # 获取每个bin的顶级ASV ID

# 过滤每个bin的主要ASVs，显著减小physeq对象的大小
physeq <- prune_taxa(bin_otu_ids, physeq)  # 保留bins的主要ASVs
# 从bins中提取主要ASVs
bins_map <- bins[, "TopTaxonID", drop = FALSE]  # 创建bins映射
# 按physeq对象中的ASVs排序
bins_map <- bins_map[match(taxa_names(physeq), bins_map$TopTaxonID), , drop = FALSE]

# 确保physeq和bins_map中的ASVs一致
identical(bins_map$TopTaxonID, taxa_names(physeq))  # 检查一致性

# 将ASVs重命名为bin名称
taxa_names(physeq) <- rownames(bins_map)  # 重命名

# 验证taxa_names或tip.labels中没有重复
length(unique(taxa_names(physeq))) == length(taxa_names(physeq))  # 检查唯一性
length(unique(phy_tree(physeq)$tip.label)) == length(phy_tree(physeq)$tip.label)  # 检查唯一性
validObject(physeq)  # 验证对象有效性
head(taxa_names(physeq))  # 查看前几个taxa名称

# 准备组装数据
bins_assembly <- as.data.frame(t(bins_assembly))  # 转置数据框
colnames(bins_assembly) <- as.character(unlist(bins_assembly[1, ]))  # 设置列名

# 转换大小写：Bin -> bin。大小写不匹配经常导致错误且容易被忽视
rownames(bins_map) <- tolower(rownames(bins_map))  # 转换为小写

# 重新排序bins_assembly以匹配bins_map中的顺序
bins_assembly <- bins_assembly[match(rownames(bins_map), rownames(bins_assembly)), , drop = FALSE]

# 检查bins_map和bins_assembly的行名是否一致
setdiff(rownames(bins_map), rownames(bins_assembly))  # 如果结果为0，则行名匹配
setdiff(rownames(bins_assembly), rownames(bins_map))  # 如果结果为0，则行名匹配

# 创建一个新的assembly数据框，将行名添加为新列
bins_assembly_new <- bins_assembly %>%
  tibble::rownames_to_column(var = "bin")  # 行名转为列
bins_assembly_new$bin <- sub("^b", "B", bins_assembly_new$bin)  # 将小写b替换为大写B

# 将列转换为数值格式
bins_assembly_new$HeS <- as.numeric(bins_assembly_new$HeS)  # 异质性选择
bins_assembly_new$HoS <- as.numeric(bins_assembly_new$HoS)  # 同质性选择
bins_assembly_new$DL <- as.numeric(bins_assembly_new$DL)    # 漂变与分散限制
bins_assembly_new$HD <- as.numeric(bins_assembly_new$HD)    # 同质性分散
bins_assembly_new$DR <- as.numeric(bins_assembly_new$DR)    # 漂变比例

# 导入bins的相对丰度数据
bins_RA <- read.csv("07Bin_Top_Taxon.csv", row.names = 1, header = TRUE)  # 读取bin的顶级分类信息
# 将行名转换为小写
rownames(bins_RA) <- tolower(rownames(bins_RA))  # 转换为小写

# 重新排序bins_RA以匹配bins_map中的顺序
bins_RA <- bins_RA[rownames(bins_map), , drop = FALSE]  # 按bins_map顺序排序
rownames(bins_RA) <- sub("^b", "B", rownames(bins_RA))  # 将小写b替换为大写B

# 检查bins_RA的行名和physeq中的taxa_names是否相同
identical(rownames(bins_RA), taxa_names(physeq))  # 检查一致性

# 创建一个新的bins_RA数据框，将行名添加为新列
bins_RA_new <- bins_RA %>%
  tibble::rownames_to_column(var = "bin")  # 行名转为列

# 筛选前10个门用于可视化准备
summarized_data <- bins_RA_new %>%
  group_by(TopTaxon.Phylum) %>%  # 按门分组
  summarise(total_RA = sum(BinRA)) %>%  # 计算每个门的总相对丰度
  ungroup() %>%
  arrange(desc(total_RA))  # 按总相对丰度降序排列

# 选择前10个门（示例数据只有两个门）
top_phyla <- summarized_data$TopTaxon.Phylum[1:10]

# 创建一个新列来指示每个bin的门是否在前10名中
bins_RA_new <- bins_RA_new %>%
  mutate(Phylum = ifelse(TopTaxon.Phylum %in% top_phyla, TopTaxon.Phylum, "others"))  # 不在前10的归为"others"

# 提取tax_table，确保它是data.frame格式
tax_table_df <- as.data.frame(physeq@tax_table)  # 转换为数据框
original_rownames <- rownames(tax_table_df)  # 保存原始行名
original_colnames <- colnames(tax_table_df)  # 保存原始列名

# 更新Phylum列，将低丰度组设置为"others"
tax_table_df$Phylum <- ifelse(tax_table_df$Phylum %in% top_phyla, tax_table_df$Phylum, "others")

# 创建一个新的physeq对象
new_physeq <- physeq  # 复制physeq对象

# 将修改后的tax_table分配回phyloseq对象
new_physeq@tax_table <- tax_table(tax_table_df)  # 更新分类表

# 恢复tax_table的行名和列名
rownames(new_physeq@tax_table) <- original_rownames  # 恢复行名
colnames(new_physeq@tax_table) <- original_colnames  # 恢复列名

# 整理bins的元数据
bins_info <- cbind(bins_assembly, bins_RA_new$BinRA)  # 合并组装信息和相对丰度
bins_info <- cbind(bins_info, bins_RA_new$Phylum)  # 添加门信息
rownames(bins_info) <- bins_assembly_new$bin  # 设置行名
colnames(bins_info) <- c("HeS", "HoS", "DL", "HD", "DR", "binsRA", "Phylum")  # 设置列名

# 从physeq对象中提取树
new_tree <- phy_tree(physeq)  # 获取系统发育树

# 将系统发育树与元数据匹配
metaplot <- match.phylo.data(new_tree, bins_info$Phylum)  # 匹配树和门信息
# 设置可视化的颜色
colors <- c("#9ACD32", "#EE6A50", "#87CEFA", "#FFC125", "#D15FEE", "#8DEEEE", "#800000",
            "#006400", "#800080", "#808080", "#B0171F", "#191970", "#7B68EE",
            "#00CD00", "Black")  # 颜色列表
phylum_colors <- setNames(colors, top_phyla)  # 将颜色映射到门

# 使用额外信息创建树
treeWithInfo <- metaplot$phy  # 获取匹配后的树

# 创建映射，将树的叶子与正确的数据关联
groupInfo1 <- split(treeWithInfo$tip.label, as.vector(bins_info$Phylum))  # 按门分组
treeWithInfo <- groupOTU(treeWithInfo, groupInfo1)  # 对树的OTU进行分组
# 确保有正确的Phylum信息
tip_data <- data.frame(
  label = treeWithInfo$tip.label,
  phylum = bins_info$Phylum[match(treeWithInfo$tip.label, rownames(bins_info))]
)
# 根据bins_map的顺序重排bins_RA_group的列
bins_RA_group <- read.csv("06Bin_relative_abundance_in_each_group_of_samples.csv", header = TRUE) 
# 从bin的列名中提取数字部分
bin_columns <- grep("^bin\\d+$", colnames(bins_RA_group), value = TRUE)
bin_numbers <- as.numeric(gsub("bin", "", bin_columns))

# 获取bins_map中bin的顺序（假设行名格式为binXX）
bins_map_numbers <- as.numeric(gsub("bin", "", rownames(bins_map)))

# 创建一个映射关系，将bins_map的顺序应用到bins_RA_group
bin_order_map <- data.frame(
  original_col = bin_columns,
  bin_number = bin_numbers,
  stringsAsFactors = FALSE
)

# 对bin_order_map按照bins_map的顺序排序
bin_order_map <- bin_order_map[order(match(bin_order_map$bin_number, bins_map_numbers)), ]

# 重新排序bins_RA_group的列
non_bin_columns <- setdiff(colnames(bins_RA_group), bin_columns)
new_column_order <- c(non_bin_columns, bin_order_map$original_col)
bins_RA_group <- bins_RA_group[, new_column_order]
# "group"变量现在是treeWithInfo中的隐藏变量
# 这已经映射到treeWithInfo中
# 当处理大型数据集时，可以注释掉`geom_tiplab`行以减少混乱
# 设置可视化的颜色
# 3.导出iTOL绘图的数据####
# 创建输出目录（如果不存在）
output_dir <- "iCAMP_visualization_data"
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

# 1. 导出系统发育树文件
library(ape)  # 用于写入树文件
# 导出最终的树文件（带有分组信息）
write.tree(treeWithInfo, file = file.path(output_dir, "iCAMP_tree.nwk"))

# 2. 导出Bin相关数据表
# 导出装配过程数据（宽格式）
write.csv(bins_assembly_new, file = file.path(output_dir, "bins_assembly_processes.csv"), row.names = FALSE)
# 导出bin相对丰度数据
write.csv(bins_RA_new, file = file.path(output_dir, "bins_relative_abundance.csv"), row.names = FALSE)
# 导出bin综合信息
write.csv(bins_info, file = file.path(output_dir, "bins_comprehensive_info.csv"))
# 将重新排序的数据保存为新文件
write.csv(bins_RA_group,file = file.path(output_dir,"bin_relative_abundance_reordered.csv") , row.names = FALSE)
# 3. 导出分类信息表
write.csv(as.data.frame(tax_table_df), file = file.path(output_dir, "taxonomy_info.csv"))
# 4.R语言可视化 ####
p1 <- ggtree(treeWithInfo,
             aes(color = group),  # 按群组着色
             layout = "fan",      # 使用扇形布局
             branch.length = 'none',  # 忽略分支长度
             open.angle = 10,     # 开放角度设置为10
             lwd = 0.5) +         # 线宽设置为0.5
  geom_tiplab(offset = 25, size = 3) +  # 添加叶子标签
  geom_tippoint(aes(color = group), size = 1) +  # 添加叶子点
  scale_color_manual(values = phylum_colors)  # 使用自定义颜色
p1  # 展示第一个图

# 将树逆时针旋转90度
p1 <- rotate_tree(p1, -90)  # 旋转树
p1  # 展示旋转后的树

# 逐层添加装配比例
p2 <- p1 + 
  geom_fruit(data = bins_assembly_new, 
             geom = geom_col,
             mapping = aes(y = bin, x = HeS), # x表示装配过程的比例
             pwidth = 0.2,  # 面板宽度
             width = 0.5,   # 条形宽度
             offset = 0.07, # 偏移量
             fill = "#BC3C29")  # 填充颜色
p2  # 展示添加了HeS后的图

p3 <- p2 + 
  geom_fruit(data = bins_assembly_new, 
             geom = geom_col,
             mapping = aes(y = bin, x = HoS),  # 添加HoS层
             pwidth = 0.2, 
             width = 0.5, 
             fill = "#0072B5")  # 填充颜色
p3  # 展示添加了HoS后的图

p4 <- p3 +
  geom_fruit(data = bins_assembly_new, 
             geom = geom_col,
             mapping = aes(y = bin, x = DL),  # 添加DL层
             pwidth = 0.2, 
             width = 0.5, 
             #grid.params=list(),
             fill = "#E18727")  # 填充颜色
p4  # 展示添加了DL后的图

p5 <- p4 + 
  geom_fruit(data = bins_assembly_new, 
             geom = geom_col,
             mapping = aes(y = bin, x = HD),  # 添加HD层
             pwidth = 0.2, 
             width = 0.5, 
             #grid.params=list(),
             fill = "#20854E")  # 填充颜色
p5  # 展示添加了HD后的图

p6 <- p5 + 
  geom_fruit(data = bins_assembly_new, 
             geom = geom_col,
             mapping = aes(y = bin, x = DR),  # 添加DR层
             pwidth = 0.2, 
             width = 0.5, 
             #grid.params=list(),
             fill = "#7876B1")  # 填充颜色
p6  # 展示添加了DR后的图

# 添加bins的相对丰度
p7 <- p6 + 
  geom_fruit(data = bins_RA_new, 
             geom = geom_col,
             mapping = aes(y = bin, x = BinRA, fill = group),  # 添加相对丰度
             pwidth = 0.5, 
             width = 0.5) +
  scale_fill_manual(values = phylum_colors)  # 使用自定义颜色
p7  # 展示最终图
# 将最终图保存为高分辨率PDF
ggsave("icamp_circle_tree.pdf", plot = p7, dpi = 1200, width = 8, height = 8)  # 保存图像











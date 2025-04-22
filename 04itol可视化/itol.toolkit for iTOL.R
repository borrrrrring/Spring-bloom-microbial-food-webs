# 加载所需的库
library(itol.toolkit)
library(data.table)
library(ape)
library(dplyr)

# 加载树文件和元数据文件
tree_1 <- read.tree("iCAMP_tree.nwk")  # 读取树文件
data_1 <-read.table("otu_table.txt",header = T) # 读取元数据

# 1.根据门（Phylum）对bin树进行着色####
unit_1 <- create_unit(
  data = data_1 %>% select(ID, Phylum),  # 选择ID和门（Phylum）列
  key = "itol_3al_2_range",  # 设置关键字
  type = "TREE_COLORS",  # 设置类型为LABEL
  color = "table2itol",
  subtype = "range",
  tree = tree_1  # 使用树文件
)
write_unit(unit_1, paste0(getwd(), "/门水平TREE_COLORS.txt"))  # 写入着色数据

# 根据门（Phylum）生成颜色条（Color Strip）
set.seed(123)  # 设置随机种子
unit_2 <- create_unit(
  data = data_1 %>% select(ID, Phylum),  # 选择ID和科（Class）列
  key = "itol_3al_3_strip",  # 设置关键字
  type = "DATASET_COLORSTRIP",  # 设置类型为颜色条
  color = "wesanderson",
  tree = tree_1  # 使用树文件
)
unit_2@common_themes$basic_theme$margin <- 50
write_unit(unit_2,paste0(getwd(),"/门水平Color Strip.txt"))  # 写入颜色条数据

# 根据不同组装过程生成多条形图（Multibar Plot）
unit_3 <- create_unit(
  data = data_1 %>% select(ID,HeS1,HoS1,DL1,HD1,DR1),  # 选择He1、2、3代表不同时期
  key = "PRE_multibar",  # 设置关键字
  type = "DATASET_MULTIBAR",  # 设置类型为多条形图
  tree = tree_1  # 使用树文件
)
unit_3@common_themes$basic_theme$margin <- 50
write_unit(unit_3, paste0(getwd(), "/PREmultibar.txt"))  # 写入多条形图数据

unit_4 <- create_unit(
  data = data_1 %>% select(ID, HeS2,HoS2,DL2,HD2,DR2),  # 选择ID、NS和OS列
  key = "IN_multibar",  # 设置关键字
  type = "DATASET_MULTIBAR",  # 设置类型为多条形图
  tree = tree_1  # 使用树文件
)
unit_4@common_themes$basic_theme$margin <- 100
write_unit(unit_4, paste0(getwd(), "/INmultibar.txt"))  # 写入多条形图数据
unit_5 <- create_unit(
  data = data_1 %>% select(ID, HeS3,HoS3,DL3,HD3,DR3),  # 选择ID、NS和OS列
  key = "POST_multibar",  # 设置关键字
  type = "DATASET_MULTIBAR",  # 设置类型为多条形图
  tree = tree_1  # 使用树文件
)
write_unit(unit_5, paste0(getwd(), "/Postmultibar.txt"))  # 写入多条形图数据


# 根据Bin丰度生成热图（Heatmap）
unit_6 <- create_unit(
  data = data_1 %>% select(ID,PREBinRA	,INEBinRA,POSTBinRA),  # 选择ID、NS和OS列
  key = "heatmap",  # 设置关键字
  type = "DATASET_HEATMAP",  # 设置类型为热图
  tree = tree_1  # 使用树文件
)
write_unit(unit_6, paste0(getwd(), "/binRAheatmap.txt"))  # 写入热图数据# 根据Bin丰度生成热图（Heatmap）
# 根据Bin丰度生成条形图
unit_7 <- create_unit(
  data = data_1 %>% select(ID,PREBinRA	,INEBinRA,POSTBinRA),  # 选择ID、NS和OS列
  key = "BinRA_multibar",  # 设置关键字
  type = "DATASET_MULTIBAR",  # 设置类型为热图
  tree = tree_1  # 使用树文件
)
write_unit(unit_7, paste0(getwd(), "/binRA_MULTIBAR.txt"))  # 写入热图数据


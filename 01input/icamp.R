rm(list=ls())   # 清除工作环境中的所有变量
t0=Sys.time()   # 记录开始时间，用于计算总耗时

# 1.数据准备####
# 存放输入文件的文件夹
wd="C:/Users/lenovo/Desktop/icamp/01input"

# OTU表文件(Tab分隔的txt文件)
com.file="otus.txt"

# 系统发育树文件
tree.file="phylo.tre"

# 分类信息(分类学)文件
clas.file="classification.txt"

# 处理信息表
treat.file="treatment.txt"

# 环境变量文件
env.file="environment.txt"
# 如果你没有环境文件或环境可能不代表生态位，跳过步骤7和8，但要检查确定分箱设置的替代方法，例如bin.size.limit

# 保存输出的文件夹。即使只是测试示例数据，也请更改为新文件夹
save.wd="C:/Users/lenovo/Desktop/icamp/output"
if(!dir.exists(save.wd)){dir.create(save.wd)}  # 如果输出文件夹不存在则创建

# 2.关键参数设置####
prefix=""   # 输出文件名的前缀，通常使用项目ID
rand.time=100   # 随机化次数，通常1000次足够。对于示例测试，可以设置为100或更少以节省时间
nworker=4       # 并行计算的线程数，取决于计算机的CPU核心数
memory.G=50     # 设置所需的内存大小(但应小于硬盘可用空间)，以便大型树的计算不受物理内存限制。单位为Gb

# 3.加载R包和数据####
library(iCAMP)  # 加载iCAMP包
library(ape)    # 加载ape包，用于处理系统发育树
setwd(wd)       # 设置工作目录为输入文件夹
comm=t(read.table(com.file, header = TRUE, sep = "\t", row.names = 1,
                  as.is = TRUE, stringsAsFactors = FALSE, comment.char = "",
                  check.names = FALSE))  # 读取OTU表并转置
tree=read.tree(file = tree.file)  # 读取系统发育树
clas=read.table(clas.file, header = TRUE, sep = "\t", row.names = 1,
                as.is = TRUE, stringsAsFactors = FALSE, comment.char = "",
                check.names = FALSE)  # 读取分类信息
treat=read.table(treat.file, header = TRUE, sep = "\t", row.names = 1,
                 as.is = TRUE, stringsAsFactors = FALSE, comment.char = "",
                 check.names = FALSE)  # 读取处理信息

env=read.table(env.file, header = TRUE, sep = "\t", row.names = 1,
               as.is = TRUE, stringsAsFactors = FALSE, comment.char = "",
               check.names = FALSE)  # 读取环境变量，如果没有env.file则跳过此步骤

# 4. 匹配OTU表和处理信息表中的样本ID####
sampid.check=match.name(rn.list=list(comm=comm,treat=treat,env=env))
# sampid.check=match.name(rn.list=list(comm=comm,treat=treat))  # 如果没有env.file，使用这行
# 对于示例数据，输出应为"All match very well"
# 对于你的数据文件，如果你还没有匹配它们的ID，不匹配的样本将被移除
treat=sampid.check$treat  # 更新处理信息
comm=sampid.check$comm    # 更新群落数据
comm=comm[,colSums(comm)>0,drop=FALSE]  # 如果某些不匹配的样本被移除，一些OTU可能成为空行，此行用于移除它们
env=sampid.check$env      # 更新环境变量，如果没有env.file则跳过此步骤

# 5.匹配OTU表和树文件中的OTU ID####
spid.check=match.name(cn.list=list(comm=comm),rn.list=list(clas=clas),tree.list=list(tree=tree))
# 对于示例数据，输出应为"All match very well"
# 对于你的数据文件，如果你之前没有匹配ID，不匹配的OTU将被移除
comm=spid.check$comm  # 更新群落数据
clas=spid.check$clas  # 更新分类信息
tree=spid.check$tree  # 更新系统发育树

# 6.计算成对系统发育距离矩阵####
# 由于微生物群落数据通常有大量的物种(OTU或ASV)，我们使用R包"bigmemory"中的"big.matrix"来处理大型系统发育距离矩阵
setwd(save.wd)  # 切换到输出文件夹
if(!file.exists("pd.desc")) 
{
  pd.big=iCAMP::pdist.big(tree = tree, wd=save.wd, nworker = nworker, memory.G = memory.G)
  # 输出文件:
  # path.rda: R对象，列出从根到每个末端的所有节点和边长。以R数据格式保存。计算系统发育距离矩阵时的中间输出
  # pd.bin: 由R包bigmemory中的big.matrix函数生成的BIN文件(backingfile)。这是存储成对系统发育距离值的大矩阵。使用这种bigmemory格式文件，在调用big matrix进行计算时，我们不需要内存而是使用硬盘
  # pd.desc: DESC文件(descriptorfile)，保存backingfile(pd.bin)的描述
  # pd.taxon.name.csv: 逗号分隔的csv文件，存储树尖端(OTU)的ID，作为大系统发育距离矩阵的行/列名
}else{
  # 如果你在之前的运行中已经计算了系统发育距离矩阵
  pd.big=list()
  pd.big$tip.label=read.csv(paste0(save.wd,"/pd.taxon.name.csv"),row.names = 1,stringsAsFactors = FALSE)[,1]
  pd.big$pd.wd=save.wd
  pd.big$pd.file="pd.desc"
  pd.big$pd.name.file="pd.taxon.name.csv"
}

#如果适用于"替代方法"(基于随机性)，你可以跳过步骤7-8
# 7.评估物种间生态位偏好差异####
# 这一步需要env变量
# 由于微生物群落数据通常有大量的物种(OTU或ASV)，我们使用R包"bigmemory"中的"big.matrix"来处理大型生态位差异矩阵
setwd(save.wd)  # 切换到输出文件夹
niche.dif=iCAMP::dniche(env = env,comm = comm,method = "niche.value",
                        nworker = nworker,out.dist=FALSE,bigmemo=TRUE,
                        nd.wd=save.wd)  # 计算生态位差异

# 8.箱内系统发育信号评估####
# 对于真实数据，你可能尝试几种不同的分箱设置，并选择导致最佳箱内系统发育信号的那一个
# 这一步需要env变量
# 8.1 # 尝试使用当前设置进行系统发育分箱
ds = 0.2  # 可以更改设置以探索最佳选择
bin.size.limit = 5  # 可以更改设置以探索最佳选择。这里设为5只是为了小型示例数据集。对于真实数据，通常尝试12到48

# taxa.binphy.big的树必须是有根树
if(!ape::is.rooted(tree))
{
  tree.rt=iCAMP::midpoint.root.big(tree = tree, pd.desc = pd.big$pd.file,
                                   pd.spname = pd.big$tip.label,pd.wd = pd.big$pd.wd,
                                   nworker = nworker)  # 如果树是无根的，使用中点法添加根
  tree=tree.rt$tree  # 更新为有根树
}
phylobin=taxa.binphy.big(tree = tree, pd.desc = pd.big$pd.file,pd.spname = pd.big$tip.label,
                         pd.wd = pd.big$pd.wd, ds = ds, bin.size.limit = bin.size.limit,
                         nworker = nworker)  # 基于系统发育距离进行分箱
# 8.2 # 测试箱内系统发育信号
sp.bin=phylobin$sp.bin[,3,drop=FALSE]  # 提取物种的箱分配
sp.ra=colMeans(comm/rowSums(comm))  # 计算物种的平均相对丰度
abcut=3  # 你可以移除一些物种，如果它们太稀少而无法执行可靠的相关测试
commc=comm[,colSums(comm)>=abcut,drop=FALSE]  # 移除稀有物种
dim(commc)  # 显示过滤后的群落矩阵维度
spname.use=colnames(commc)  # 获取要使用的物种名称
binps=iCAMP::ps.bin(sp.bin = sp.bin,sp.ra = sp.ra,spname.use = spname.use,
                    pd.desc = pd.big$pd.file, pd.spname = pd.big$tip.label, pd.wd = pd.big$pd.wd,
                    nd.list = niche.dif$nd,nd.spname = niche.dif$names,ndbig.wd = niche.dif$nd.wd,
                    cor.method = "pearson",r.cut = 0.1, p.cut = 0.05, min.spn = 5)  # 计算箱内系统发育信号
if(file.exists(paste0(prefix,".PhyloSignalSummary.csv"))){appendy=TRUE;col.namesy=FALSE}else{appendy=FALSE;col.namesy=TRUE}
write.table(data.frame(ds=ds,n.min=bin.size.limit,binps$Index),file = paste0(prefix,".PhyloSignalSummary.csv"),
            append = appendy, quote=FALSE, sep=",", row.names = FALSE,col.names = col.namesy)  # 将系统发育信号摘要写入CSV文件
if(file.exists(paste0(prefix,".PhyloSignalDetail.csv"))){appendy2=TRUE;col.namesy2=FALSE}else{appendy2=FALSE;col.namesy2=TRUE}
write.table(data.frame(ds=ds,n.min=bin.size.limit,binID=rownames(binps$detail),binps$detail),file = paste0(prefix,".PhyloSignalDetail.csv"),
            append = appendy2, quote = FALSE, sep = ",", row.names = FALSE, col.names = col.namesy2)  # 将系统发育信号详情写入CSV文件
# 由于这个示例小数据是随机生成的，相关性应该非常弱
# 通常，你需要寻找一个能导致更高RAsig.abj(具有显著系统发育信号的箱的相对丰度)和相对高meanR(箱间的平均相关系数)的分箱设置
# 参见函数"ps.bin"的帮助文档，了解输出的含义

# 9. iCAMP分析####
# 9.1 # 不省略小箱
# 常用 # 将sig.index设置为Confidence而不是SES.RC (betaNRI/NTI + RCbray)
bin.size.limit = 24 # 对于真实数据，通常根据系统发育信号测试使用适当的数字或尝试一些设置然后选择合理的随机性水平。我们的经验是12、24或48。但对于这个太小的示例数据集，必须使用5
sig.index="Confidence"  # 在icamp.big的帮助文档中查看其他选项
icres=iCAMP::icamp.big(comm=comm, pd.desc = pd.big$pd.file, pd.spname=pd.big$tip.label,
                       pd.wd = pd.big$pd.wd, rand = rand.time, tree=tree,
                       prefix = prefix, ds = 0.2, pd.cut = NA, sp.check = TRUE,
                       phylo.rand.scale = "within.bin", taxa.rand.scale = "across.all",
                       phylo.metric = "bMPD", sig.index=sig.index, bin.size.limit = bin.size.limit, 
                       nworker = nworker, memory.G = memory.G, rtree.save = FALSE, detail.save = TRUE, 
                       qp.save = FALSE, detail.null = FALSE, ignore.zero = TRUE, output.wd = save.wd, 
                       correct.special = TRUE, unit.sum = rowSums(comm), special.method = "depend",
                       ses.cut = 1.96, rc.cut = 0.95, conf.cut=0.975, omit.option = "no",meta.ab = NULL)  # 进行iCAMP分析
# 该函数有很多参数，请查看"icamp.big"的帮助文档
# 输出文件:
# Test.iCAMP.detail.rda: 以R数据格式保存的对象"icres"。它是一个列表对象。第一个元素bNRIiRCa是每个配对比较(每个周转)中每个组装过程的相对重要性结果。第二个元素"detail"包括分箱信息(命名为taxabin)、每个箱中的系统发育和分类度量结果(命名为bNRIi、RCa等)、每个箱的相对丰度(bin.weight)、群落间每个周转中每个过程的相对重要性(processes)、输入设置(setting)和输入群落数据矩阵(comm)。有关更多详细信息，请参见函数icamp.big的帮助文档

# 9.2至9.4是一些可选的特殊设置，你可以探索
# 9.2 # 探索不同的零模型显著性测试方法
# 9.2.1 # 设置detail.null=TRUE，输出所有零值，以便进行正态性测试并在不同选项之间切换
detail.null=TRUE
bin.size.limit = 5 
sig.index="SES.RC"  # 这是传统方法，假设系统发育度量的零值遵循正态分布
prefixb="TestB"

icres2=iCAMP::icamp.big(comm=comm, pd.desc = pd.big$pd.file, pd.spname=pd.big$tip.label,
                        pd.wd = pd.big$pd.wd, rand = rand.time, tree=tree,
                        prefix = prefixb, ds = 0.2, pd.cut = NA, sp.check = TRUE,
                        phylo.rand.scale = "within.bin", taxa.rand.scale = "across.all",
                        phylo.metric = "bMPD", sig.index=sig.index, bin.size.limit = bin.size.limit, 
                        nworker = nworker, memory.G = memory.G, rtree.save = FALSE, detail.save = TRUE, 
                        qp.save = FALSE, detail.null = detail.null, ignore.zero = TRUE, output.wd = save.wd, 
                        correct.special = TRUE, unit.sum = rowSums(comm), special.method = "depend",
                        ses.cut = 1.96, rc.cut = 0.95, conf.cut=0.975, omit.option = "no",meta.ab = NULL)  # 使用SES.RC进行iCAMP分析
# 9.2.2 # 正态性测试
nntest=iCAMP::null.norm(icamp.output=icres2, p.norm.cut=0.05, detail.out=FALSE)  # 对零模型值进行正态性测试
# 输出显示每个箱中非正态分布的比例，即具有显著偏离正态分布的零值的周转比例
# 如果某些比例值很高，可能需要更改为使用"Confidence"作为sig.index

# 9.2.3 # 将sig.index更改为"Confidence"
icres3=iCAMP::change.sigindex(icamp.output = icres2, sig.index = "Confidence", detail.save = TRUE, detail.null = FALSE, conf.cut = 0.975)  # 将显著性指数更改为Confidence
head(icres3$CbMPDiCBraya)  # 显示结果的前几行

# 9.2.4 # 将sig.index更改为"RC"，同时用于系统发育和分类度量
icres4=iCAMP::change.sigindex(icamp.output = icres2, sig.index = "RC", detail.save = TRUE, detail.null = FALSE, rc.cut = 0.95)  # 将显著性指数更改为RC
head(icres4$RCbMPDiRCbraya)  # 显示结果的前几行

# 9.2.5 # 该函数还可以更改显著性阈值
icres5=iCAMP::change.sigindex(icamp.output = icres2, sig.index = "SES.RC", detail.save = TRUE, detail.null = FALSE, ses.cut = 1.64, rc.cut = 0.9)  # 更改显著性阈值
head(icres5$bNRIiRCbraya)  # 显示结果的前几行

# 9.3 # 如果你输入的"comm"中的平均相对丰度与区域库中不同，你可以指定区域库中每个物种的相对丰度
meta.ab=rep(1,ncol(comm))  # 这里假设区域库中所有物种实际上具有相同的相对丰度
prefix2=paste0(prefix,".MetaCrct")
sig.index="Confidence"  # 查看icamp.big的帮助文档中的其他选项
icres.meta=iCAMP::icamp.big(comm=comm, pd.desc = pd.big$pd.file, pd.spname=pd.big$tip.label,
                            pd.wd = pd.big$pd.wd, rand = rand.time, tree=tree,
                            prefix = prefix2, ds = 0.2, pd.cut = NA, sp.check = TRUE,
                            phylo.rand.scale = "within.bin", taxa.rand.scale = "across.all",
                            phylo.metric = "bMPD", sig.index=sig.index, bin.size.limit = bin.size.limit, 
                            nworker = nworker, memory.G = memory.G, rtree.save = FALSE, detail.save = TRUE, 
                            qp.save = FALSE, detail.null = FALSE, ignore.zero = TRUE, output.wd = save.wd, 
                            correct.special = TRUE, unit.sum = rowSums(comm), special.method = "depend",
                            ses.cut = 1.96, rc.cut = 0.95, conf.cut=0.975, omit.option = "no",meta.ab=meta.ab)  # 使用指定的元群落丰度进行iCAMP分析

# 9.4 # 考虑省略小箱
# 9.4.1 # 如果你想省略小箱而不是将它们合并到最近的亲缘中，将omit.option设置为"test"以检查将省略哪些内容
omit.option = "test"
icres.omit=iCAMP::icamp.big(comm=comm, pd.desc = pd.big$pd.file, pd.spname=pd.big$tip.label,
                            pd.wd = pd.big$pd.wd, rand = rand.time, tree=tree,
                            prefix = prefix, ds = 0.2, pd.cut = NA, sp.check = TRUE,
                            phylo.rand.scale = "within.bin", taxa.rand.scale = "across.all",
                            phylo.metric = "bMPD", sig.index=sig.index, bin.size.limit = bin.size.limit, 
                            nworker = nworker, memory.G = memory.G, rtree.save = FALSE, detail.save = TRUE, 
                            qp.save = FALSE, detail.null=FALSE, ignore.zero = TRUE, output.wd = save.wd, 
                            correct.special = TRUE, unit.sum = rowSums(comm), special.method = "depend",
                            ses.cut = 1.96, rc.cut = 0.95, conf.cut=0.975, omit.option = omit.option)  # 测试哪些箱将被省略
# "test"将返回被省略物种的详细表格

# 9.4.2 # 然后将其设置为"omit"以省略小箱
omit.option = "omit"
icres.omit2=iCAMP::icamp.big(comm=comm, pd.desc = pd.big$pd.file, pd.spname=pd.big$tip.label,
                             pd.wd = pd.big$pd.wd, rand = rand.time, tree=tree,
                             prefix = prefix, ds = 0.2, pd.cut = NA, sp.check = TRUE,
                             phylo.rand.scale = "within.bin", taxa.rand.scale = "across.all",
                             phylo.metric = "bMPD", sig.index=sig.index, bin.size.limit = bin.size.limit, 
                             nworker = nworker, memory.G = memory.G, rtree.save = FALSE, detail.save = TRUE, 
                             qp.save = FALSE, detail.null=FALSE, ignore.zero = TRUE, output.wd = save.wd, 
                             correct.special = TRUE, unit.sum = rowSums(comm), special.method = "depend",
                             ses.cut = 1.96, rc.cut = 0.95, conf.cut=0.975, omit.option = omit.option)  # 省略小箱
# 在这个简单示例中，由于所有箱都很小，"omit"应该返回错误。在真实数据中，这将继续对足够大的严格箱(>bin.size.limit)进行iCAMP分析

# 9.5 # 输入相对丰度的群落矩阵(值<1)而不是计数
comra=comm/rowSums(comm)  # 将OTU计数转换为相对丰度
prefixra=paste0(prefix,"RA")
bin.size.limit = 5  # 对于真实数据，通常根据系统发育信号测试使用适当的数字或尝试一些设置然后选择合理的随机性水平。我们的经验是12、24或48。但对于这个太小的示例数据集，必须使用5
icres6=iCAMP::icamp.big(comm=comra,tree=tree,pd.desc=pd.big$pd.file, pd.spname=pd.big$tip.label, pd.wd=pd.big$pd.wd,
                        rand=rand.time,prefix=prefixra,ds=0.2,pd.cut=NA,sp.check=TRUE,
                        phylo.rand.scale="within.bin",taxa.rand.scale="across.all",
                        phylo.metric="bMPD",sig.index="Confidence",
                        bin.size.limit=bin.size.limit,nworker=nworker,memory.G=memory.G,
                        rtree.save=FALSE,detail.save=TRUE,qp.save=FALSE,detail.null=FALSE,
                        ignore.zero=TRUE,output.wd=save.wd,correct.special=TRUE,unit.sum=rowSums(comra),
                        special.method="depend",ses.cut = 1.96,rc.cut = 0.95,conf.cut=0.975,
                        omit.option="no",meta.ab=NULL, taxo.metric="bray", transform.method=NULL,
                        logbase=2, dirichlet=TRUE)  # 使用相对丰度数据进行iCAMP分析

# 9.6 # 群落数据转换和分类相异性指数更改
taxo.metric='euclidean'  # 使用欧式距离作为分类相异性指数
transform.method='hellinger'  # 使用Hellinger转换
prefixtran=paste0(prefix,"Hel")
bin.size.limit = 24  # 对于真实数据，通常根据系统发育信号测试使用适当的数字或尝试一些设置然后选择合理的随机性水平。我们的经验是12、24或48。但对于这个太小的示例数据集，必须使用5
icres7=iCAMP::icamp.big(comm=comm,tree=tree,pd.desc=pd.big$pd.file, pd.spname=pd.big$tip.label, pd.wd=pd.big$pd.wd,
                        rand=rand.time,prefix=prefixtran,ds=0.2,pd.cut=NA,sp.check=TRUE,
                        phylo.rand.scale="within.bin",taxa.rand.scale="across.all",
                        phylo.metric="bMPD",sig.index="Confidence",
                        bin.size.limit=bin.size.limit,nworker=nworker,memory.G=memory.G,
                        rtree.save=FALSE,detail.save=TRUE,qp.save=FALSE,detail.null=FALSE,
                        ignore.zero=TRUE,output.wd=save.wd,correct.special=TRUE,unit.sum=rowSums(comra),
                        special.method="depend",ses.cut = 1.96,rc.cut = 0.95,conf.cut=0.975,
                        omit.option="no",meta.ab=NULL, taxo.metric=taxo.metric, transform.method=transform.method,
                        logbase=2, dirichlet=FALSE)  # 使用Hellinger转换和欧式距离进行iCAMP分析

# 10.iCAMP bin级别统计####
# 对iCAMP结果按照物种组(bin)进行汇总分析
icbin=iCAMP::icamp.bins(icamp.detail = icres$detail,treat = treat,
                        clas=clas,silent=FALSE, boot = TRUE,
                        rand.time = rand.time,between.group = TRUE)
save(icbin,file = paste0(prefix,".iCAMP.Summary.rda")) # 将结果保存为R数据文件，会自动压缩并便于后续导入R
# 以下是将多个结果输出为CSV文件
write.csv(icbin$Ptuv,file = paste0("01ProcessImportance_EachTurnover.csv"),row.names = FALSE) # 每对样本间周转中各生态过程的相对重要性
write.csv(icbin$Ptk,file = paste0("02ProcessImportance_EachBin_EachGroup.csv"),row.names = FALSE) # 每个样本组中每个bin的各生态过程相对重要性
write.csv(icbin$Pt,file = paste0("03ProcessImportance_EachGroup.csv"),row.names = FALSE) # 每个样本组中各生态过程的相对重要性
write.csv(icbin$BPtk,file = paste0("04BinContributeToProcess_EachGroup.csv"),row.names = FALSE) # 每个bin对每个样本组中各生态过程的贡献
write.csv(icbin$BRPtk, file = paste0("05Bin_relative_contribution_to_each_process.csv"), row.names = FALSE) # bin对不同过程的相对贡献
write.csv(icbin$Binwt, file = paste0("06Bin_relative_abundance_in_each_group_of_samples.csv"), row.names = FALSE)#bin在不同组的相对丰度
write.csv(icbin$Bin.TopClass,file = paste0("07Bin_Top_Taxon.csv"),row.names = FALSE) # 每个bin的丰度、主要物种及其分类信息
write.csv(icbin$Class.Bin, file = paste0("08Bin_ID_classification_information.csv"), row.names = FALSE)
write.csv(data.frame(ID=rownames(icbin$Class.Bin),icbin$Class.Bin,stringsAsFactors = FALSE),
          file = paste0("Taxon_Bin.csv"),row.names = FALSE) # 每个物种的bin ID和分类信息
# 输出文件说明:详情?icamp.bins
# Test.iCAMP.Summary.rda：以R数据格式保存的"icbin"对象，详细内容可参见icamp.bins函数的帮助文档
# Test.ProcessImportance_EachGroup.csv：各生态过程在群落周转中的相对重要性
# Test.ProcessImportance_EachBin_EachGroup.csv：各生态过程在每个bin的群落周转中的相对重要性
# Test.ProcessImportance_EachTurnover.csv：各生态过程在每对样本间周转中的相对重要性
# Test.BinContributeToProcess_EachGroup.csv：每个bin对各生态过程的贡献度
# Test.Taxon_Bin.csv：显示每个物种的bin ID和分类信息的矩阵
# Test.Bin_TopTaxon.csv：显示每个bin的相对丰度、主要物种ID、在bin中的百分比及其分类信息

# 11.自举检验(Bootstrapping test)####
# 指定处理信息表中的哪一列进行分析
i=1
treat.use=treat[,i,drop=FALSE]
icamp.result=icres$CbMPDiCBraya
# 进行自举检验，评估不同组间生态过程差异的显著性
icboot=iCAMP::icamp.boot(icamp.result = icamp.result,treat = treat.use,rand.time = rand.time,
                         compare = TRUE,silent = FALSE,between.group = TRUE,ST.estimation = TRUE)
save(icboot,file=paste0(prefix,".iCAMP.Boot.",colnames(treat)[i],".rda")) # 保存自举检验结果
write.csv(icboot$summary,file = paste0(prefix,".iCAMP.BootSummary.",colnames(treat)[i],".csv"),row.names = FALSE) # 保存自举检验摘要
write.csv(icboot$compare,file = paste0(prefix,".iCAMP.Compare.",colnames(treat)[i],".csv"),row.names = FALSE) # 保存组间比较结果

# 输出文件说明:
# Test.iCAMP.Boot.Management.rda：以R数据格式保存的"icboot"对象
# Test.BootSummary.Management.csv：自举检验结果的摘要表
# Test.Compare.Management.csv：各组间比较指数、效应量和显著性水平的摘要表

# 12.其他方法: QPEN (基于整个群落的零模型分析来量化群落构建过程)####
# 12.1 # QPEN计算
qpout=iCAMP::qpen(comm=comm,pd=pd.big$pd.file,pd.big.wd=pd.big$pd.wd,
                  pd.big.spname=pd.big$tip.label,ab.weight=TRUE,
                  rand.time=rand.time, nworker=nworker,project=prefix,
                  wd=save.wd, save.bNTIRC=TRUE)
# 12.2 # 显著性检验
qptest=qpen.test(qpen.result = qpout,treat = treat,rand.time = rand.time,
                 between.group = TRUE,out.detail=TRUE,silent=FALSE)
write.csv(qptest$obs.summary,file = paste0(prefix,".QPEN.Index.Obs.Summary.csv"),row.names = FALSE) # 保存观测指数摘要
write.csv(qptest$boot.summary,file = paste0(prefix,".QPEN.Bootstrapping.Summary.csv"),row.names = FALSE) # 保存自举检验摘要
write.csv(qptest$compare,file = paste0(prefix,".QPEN.Comparison.Summary.csv"),row.names = FALSE) # 保存比较摘要
save(qptest,file = paste0(prefix,".QPEN.bootstrap.rda")) # 保存完整的检验结果

# 13.其他方法: 中性物种百分比分析####
# 分析符合中性理论的物种比例
snmout=iCAMP::snm.comm(comm = comm, treat = treat, 
                       rand=rand.time, alpha=0.05)
write.csv(snmout$stats,file = paste0(prefix,".NeutralModel.Stats.csv")) # 保存中性模型统计结果
write.csv(snmout$ratio.summary,file = paste0(prefix,".NeutralModel.TypeRatio.csv")) # 保存各类型物种比例摘要

# 14.其他方法: tNST和pNST (分类学和系统发育归一化随机性比率)####
# 如果尚未安装NST包，则安装该包
if(!("NST" %in% installed.packages()[,"Package"])){install.packages("NST")}
library(NST)
i=1
treat.use=treat[,i,drop=FALSE]

# 14.1a # 计算tNST (分类学归一化随机性比率)
tnstout=NST::tNST(comm=comm, group=treat.use, dist.method="bray", 
                  abundance.weighted=TRUE, rand=rand.time,  
                  nworker=nworker, null.model="PF", output.rand = TRUE,
                  SES = TRUE, RC = TRUE)
write.csv(tnstout$index.grp,file = paste0(prefix,".tNST.summary.",colnames(treat)[i],".csv")) # 保存tNST摘要结果
write.csv(tnstout$index.pair.grp,file = paste0(prefix,".tNST.pairwise.",colnames(treat)[i],".csv")) # 保存成对tNST结果

# 14.1b # tNST的自举检验
tnst.bt=NST::nst.boot(nst.result=tnstout, group=treat.use,
                      rand=rand.time, nworker=nworker)
write.csv(tnst.bt$NST.summary,file = paste0(prefix,".tNST.bootstr.",colnames(treat)[i],".csv")) # 保存tNST自举检验摘要
write.csv(tnst.bt$NST.compare,file = paste0(prefix,".tNST.compare.",colnames(treat)[i],".csv")) # 保存tNST组间比较结果

# 14.2a # 计算pNST (系统发育归一化随机性比率)
pnstout=NST::pNST(comm=comm, pd.desc=pd.big$pd.file, pd.wd=pd.big$pd.wd, 
                  pd.spname=pd.big$tip.label, group=treat.use, abundance.weighted=TRUE,
                  rand=rand.time, phylo.shuffle=TRUE, nworker=nworker,
                  output.rand = TRUE, SES=FALSE, RC=FALSE)
write.csv(pnstout$index.grp,file = paste0(prefix,".pNST.summary.",colnames(treat)[i],".csv")) # 保存pNST摘要结果
write.csv(pnstout$index.pair.grp,file = paste0(prefix,".pNST.pairwise.",colnames(treat)[i],".csv")) # 保存成对pNST结果

# 14.2b # pNST的自举检验
pnst.bt=NST::nst.boot(nst.result=pnstout, group=treat.use,
                      rand=rand.time, nworker=nworker)
write.csv(pnst.bt$NST.summary,file = paste0(prefix,".pNST.bootstr.",colnames(treat)[i],".csv")) # 保存pNST自举检验摘要
write.csv(pnst.bt$NST.compare,file = paste0(prefix,".pNST.compare.",colnames(treat)[i],".csv")) # 保存pNST组间比较结果

# 15.汇总核心物种、稀有物种和其他物种####
# 15.1 # 从category.txt文件中定义不同物种的类型
setwd(wd)
cate.file="category.txt"
cate=read.table(cate.file, header = TRUE, sep = "\t", row.names = 1,
                as.is = TRUE, stringsAsFactors = FALSE, comment.char = "",
                check.names = FALSE)
cate=cate[which(rownames(cate) %in% colnames(comm)),,drop=FALSE] # 移除不匹配的物种
setwd(save.wd)

# 15.2 # 按照不同物种类别分析iCAMP结果
iccate=icamp.cate(icamp.bins.result = icbin,comm = comm,cate = cate,
                  treat = treat, silent = FALSE,between.group = TRUE)
write.csv(iccate$Ptuvx,file = paste0(prefix,".iCAMP.Process_EachTurnover_EachCategory.csv")) # 保存每对样本中每类物种的生态过程重要性
write.csv(iccate$Ptx,file = paste0(prefix,".iCAMP.Process_EachGroup_EachCategory.csv")) # 保存每个样本组中每类物种的生态过程重要性

(t=format(Sys.time()-t0)) # 计算代码执行所需的时间
# 结束 #
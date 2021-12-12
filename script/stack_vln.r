library(Seurat)
library(ggplot2)
library(pathcwork)
library(tidyverse)


gene=read.table("gene.txt",header=T,sep="\t",quote='',stringsAsFactors = F)
obj=readRDS('obj_5000.rds')
top1=gene%>%group_by(group)%>%top_n(n=3,wt=y)

#从Seurat对象中提取细胞注释以及基因表达量
vln.dat=FetchData(obj,c(top1$gene,"Major.CellType","Minor.CellType"))

# 定义因子顺序，防止画图时对细胞注释进行重排
vln.dat$Major.CellType=factor(vln.dat$Major.CellType,levels = c("T/NK cells","B lymphocytes",
                                                                          "Myeloid cells","MAST cells",
                                                                          "Endothelial cells","Epithelial cells","Fibroblasts"))
vln.dat=vln.dat[order(vln.dat$Major.CellType),]
vln.dat$Minor.CellType=factor(vln.dat$Minor.CellType,levels = unique(vln.dat$Minor.CellType))
vln.dat.melt=vln.dat%>%reshape2::melt(,top1$gene)%>%rename("Gene"="variable")%>%group_by(Minor.CellType,Gene)%>%mutate(fillcolor=mean(value))

# 小提琴图的填充颜色
pal=colorRampPalette(RColorBrewer::brewer.pal(n = 9,name="YlOrRd"))(100)

# 堆积小提琴图
p1 = ggplot(vln.dat.melt,aes(x=Minor.CellType,y=value,fill=fillcolor))+
  # 把小提琴图的外缘轮廓去除
  geom_violin(linetype="blank",scale = "width")+
  scale_fill_gradientn(colors=pal,name="Average Expression")+
  # 分面
  facet_grid(Gene~.)+
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.text.x = element_text(angle = 90,hjust = 1,vjust = .5),
        strip.text.y = element_text(angle = 0),
        axis.title.x = element_blank(),
        axis.ticks = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        strip.background = element_blank(),
        legend.position = "left"
        
        )
# 画最右边的基因注释线段，这里我们用geom_segment来画，当然也可以用geom_tile()和geom_bar()来画，大家可以自己探索一下
a=top1%>%count(group)
# 用geom_segment，我们需要计算线段的起始点和终止点
# 计算每个基因group的成员数量，然后进行累加，每个累加值就是线段的终点，最后减去每个组的成员数量，就是线段的起始位点
a$yend=cumsum(a$n)
a$ystart=a$yend-a$n
# 计算y轴的text的位置，就是group名称的位置，应该是每个线段的中点
a$label_position=a$ystart+(a$yend-a$ystart)/2
p2=ggplot(a,aes(x=1,y=ystart,color=group))+
  geom_segment(aes(xend=1,yend=yend),size=3)+
  scale_color_brewer(palette = "Set1")+
  scale_y_continuous(position = "right",
                     breaks = a$label_position,
                     labels = a$group,expand = c(0,0))+
  scale_x_continuous(expand = c(0,0))+
  facet_grid(group~.,scales = 'free')+
  theme(strip.text = element_blank(),
        strip.background = element_blank(),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "none")

# 我们用geom_tile()完成细胞注释，当然也可以向上面一样用geom_segment或者geom_bar来完成
p3=ggplot(vln.dat%>%
         select(Major.CellType,Minor.CellType)%>%
         unique()%>%
         mutate(value='A'),
       aes(x=Minor.CellType,y=value,fill=Major.CellType))+
  geom_tile()+
  scale_fill_brewer(palette = "Set2")+
  scale_x_discrete(expand = c(0,0))+
  scale_y_discrete(expand = c(0,0))+
  facet_grid(.~Major.CellType,scales = "free",space = "free",switch = 'x')+
  theme(panel.background = element_blank(),
        strip.text = element_text(angle = 90,hjust = 1,vjust = .5),
        strip.background = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        legend.position = "none")

p1+p2+p3+plot_spacer()+plot_layout(ncol = 2, widths  = c(3, .03),heights = c(3,.03))
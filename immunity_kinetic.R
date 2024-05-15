library(deSolve)
library(ggplot2)
library(gtable)
library(grid)


hinvert_title_grob <- function(grob){
  # 交换宽度
  widths <- grob$widths
  grob$widths[1] <- widths[3]
  grob$widths[3] <- widths[1]
  grob$vp[[1]]$layout$widths[1] <- widths[3]
  grob$vp[[1]]$layout$widths[3] <- widths[1]
  
  # 修改对齐
  grob$children[[1]]$hjust <- 1 - grob$children[[1]]$hjust 
  grob$children[[1]]$vjust <- 1 - grob$children[[1]]$vjust 
  grob$children[[1]]$x <- unit(1, "npc") - grob$children[[1]]$x
  grob
}

add_yaxis_left <- function(g1, g2) {
  # 将坐标轴添加到左侧
  # 添加坐标轴
  pos <- c(subset(g1$layout, name == "ylab-l", select = t:r))
  index <- which(g2$layout$name == "axis-l")
  yaxis <- g2$grobs[[index]]
  # 先添加 3mm 间距
  g <- gtable_add_cols(g1, unit(3, "mm"), pos$l - 1)
  # 再添加轴
  g <- gtable_add_cols(g, g2$widths[g2$layout[index, ]$l], pos$l - 1)
  g <- gtable_add_grob(g, yaxis, pos$t, pos$l, pos$b, pos$l, clip = "off", name = "axis-l")
  # 添加轴标签
  index <- which(g2$layout$name == "ylab-l")
  ylab <- g2$grobs[[index]]
  g <- gtable_add_cols(g, g2$widths[g2$layout[index, ]$l], pos$l - 1)
  g <- gtable_add_grob(g, ylab, pos$t, pos$l, pos$b, pos$l, clip = "off", name = "ylab-l")
  g
}

add_yaxis_right <- function(g1, g2, pos) {
  # 将坐标轴添加到右侧
  # ============ 2. 轴标签 ============ #
  index <- which(g2$layout$name == "ylab-l")
  ylab <- g2$grobs[[index]]
  ylab <- hinvert_title_grob(ylab)
  # 添加轴标签
  g <- gtable_add_cols(g1, g2$widths[g2$layout[index, ]$l], pos$r)
  g <- gtable_add_grob(g, ylab, pos$t, pos$r + 1, pos$b, pos$r + 1, clip = "off", name = "ylab-r")
  # ============ 3. 轴设置 ============ #
  index <- which(g2$layout$name == "axis-l")
  yaxis <- g2$grobs[[index]]
  # 将 Y 轴线移动到最左边
  yaxis$children[[1]]$x <- unit.c(unit(0, "npc"), unit(0, "npc"))
  # 交换刻度线和刻度标签
  ticks <- yaxis$children[[2]]
  ticks$widths <- rev(ticks$widths)
  ticks$grobs <- rev(ticks$grobs)
  # 移动刻度线
  ticks$grobs[[1]]$x <- ticks$grobs[[1]]$x - unit(1, "npc") + unit(3, "pt")
  # 刻度标签位置转换和对齐
  ticks$grobs[[2]] <- hinvert_title_grob(ticks$grobs[[2]])
  yaxis$children[[2]] <- ticks
  # 添加轴，unit(3, "mm") 增加轴间距
  g <- gtable_add_cols(g, g2$widths[g2$layout[index, ]$l] + unit(3, "mm"), pos$r)
  g <- gtable_add_grob(g, yaxis, pos$t, pos$r + 1, pos$b, pos$r + 1, clip = "off", name = "axis-r")
  g
}

add_yaxis <- function(g1, g2, offset = 0) {
  # ============ 1. 主绘图区 ============ #
  # 获取主绘图区域
  pos <- c(subset(g1$layout, name == "panel", select = t:r))
  # 添加图形
  g <- gtable_add_grob(g1, g2$grobs[[which(g2$layout$name == "panel")]], 
                       pos$t, pos$l, pos$b * ((offset - 2) * 0.00001 + 1), pos$l)
  if (offset > 3 && offset %% 2 == 0) {
    g <- add_yaxis_left(g, g2)
  } else {
    g <- add_yaxis_right(g, g2, pos)
  }
  g
}

# 接受可变参数，可添加多个 Y 轴
plot_multi_yaxis <- function(..., right_label_reverse = TRUE) {
  args <- list(...)
  my_theme <- theme(panel.grid = element_blank(), panel.background = element_rect(fill = NA))
  len <- length(args)
  args[[1]] <- args[[1]] + my_theme
  g <- ggplotGrob(args[[1]])
  for (i in len:2) { 
    if (i < 4 || i %% 2 && right_label_reverse) {
      # 为轴标签添加旋转
      args[[i]] <- args[[i]] + 
        theme(axis.title.y = element_text(angle = 270))
    }
    args[[i]] <- args[[i]] + my_theme
    # 获取 gtable 对象
    g2 <- ggplotGrob(args[[i]])
    g <- add_yaxis(g, g2, offset = i)
  }
  # 绘制图形
  grid.newpage()
  grid.draw(g)
}

####################################### Immunity non-existance

r = 1
phi = 5*10^-8
eps = 0
KD = 2.2*10^6
KC = 10^9

beta = 100
phi = 5*10^-8
omg = 1

alpha = 0.97
KI = 2.4*10^7
KN = 10^5
Iomg = 0.05



derivs <- function(t,state, parms) {
  with(as.list(c(state, parms)), {
    
    dB <- r*B*(1-B/KC) - B*phi*P - eps*I*B/(1+B/KD)
    
    dP <- beta*phi*P*B - omg*P
    
    dI <- alpha * I * (1-I/KI) * B/(B + KN) - Iomg*I
    
    list(c(dB,dP,dI))
  })
}

state <- c(B = 10^6,P = 10^7,I = 0)
times <- seq(0, 25, 0.1)

parms <- c(r,phi,eps,KD,KC,
           beta,phi,omg,
           alpha,KI,KN)

yout <- ode(y = state, times = times, func = derivs, parms = parms,method = 'ode45')
#plot(yout)

########################## visualization
library(ggplot2)
options(scipen = 200)

dat = data.frame(Time = rep(yout[,1],3),
                 Value = c(yout[,2],yout[,3],yout[,4]),
                 Type = c(rep('Bacteria',nrow(yout)),
                          rep('Phage',nrow(yout)),
                          rep('Immunity cells',nrow(yout))
                 )
)


dat_b = subset(dat,dat$Type=='Bacteria')
dat_p = subset(dat,dat$Type=='Phage')
dat_i = subset(dat,dat$Type=='Immunity cells')

library(ggplot2)

p1 = ggplot(data=dat_b,aes(x=Time,y=Value))+
  geom_line(size=1.5,color="#FF6600") + theme_classic() + 
  theme(legend.position = 'none') + 
  ylab('CFU/g') + ggtitle('Bacteria') + xlab('') +
  scale_y_continuous(labels = c(expression(italic(0)),
                                expression(2.5%*%10^5),
                                expression(5.0%*%10^5),
                                expression(7.5%*%10^5),
                                expression(10.0%*%10^5)),
                     position = "left",
                     expand = c(0,0),
                     breaks = c(0,250000,500000,750000,1000000),
                     limits = c(0,1100000))


p2 = ggplot(data=dat_p,aes(x=Time,y=Value))+
  geom_line(size=1.5,color="#33CC99") + theme_classic() + 
  theme(legend.position = 'none') + 
  ylab('PFU/g') + ggtitle('Phage') + xlab('') +
  scale_y_continuous(labels = c(expression(italic(0)),
                                expression(2.5%*%10^7),
                                expression(5.0%*%10^7),
                                expression(7.5%*%10^7),
                                expression(10.0%*%10^7)),
                     position = "left",
                     expand = c(0,0),
                     breaks = c(0,25000000,50000000,75000000,100000000),
                     limits = c(0,110000000))


p3 = ggplot(data=dat_i,aes(x=Time,y=Value))+
  geom_line(size=1.5,color="#6633FF") + theme_classic() + 
  theme(legend.position = 'none') + 
  ylab('Cells/g') + ggtitle('Immunity cells') + xlab('')



dd = data.frame(cbind(Time = yout[,1],Bacteria = yout[,2],Phage = yout[,3],Immunity = yout[,4]))

colors <- c('#E83529', '#106EB8', '#F39800')
p1 <- ggplot(dd, aes(Time, Bacteria, group = 1)) + 
  geom_line(size=1.5,colour = colors[1], position = position_nudge(x = -0.2)) + 
  labs(x = "Time(h)", y = "Bacteria(CFU/g)") +
  scale_y_continuous(limits = c(0,1200000), expand = c(0,0),
                     labels = c(expression(italic(0)),
                                expression(2.5%*%10^5),
                                expression(5.0%*%10^5),
                                expression(7.5%*%10^5),
                                expression(10.0%*%10^5)),
                     breaks = c(0,250000,500000,750000,1000000)) +
  theme(panel.grid = element_blank(), 
        panel.background = element_rect(fill = NA), 
        axis.text.y = element_text(color = colors[1],size = 20),
        axis.ticks.y = element_line(color = colors[1]),
        axis.title.x = element_text(size = 23),
        axis.title.y = element_text(color = colors[1],size = 20), 
        axis.line.y = element_line(color = colors[1]),
        axis.line.x = element_line(color = 'black'),
        axis.text.x = element_text(size = 20,angle = 45, hjust = 1, vjust = 1)
  )


p2 <- ggplot(dd, aes(Time, Phage, group = 1)) + 
  geom_line(size=1.5,colour = colors[2]) + 
  labs(x = "Time(h)", y = "Phage(PFU/g)") +
  scale_y_continuous(limits = c(0,120000000), expand = c(0,0),
                     labels = c(expression(italic(0)),
                                expression(2.5%*%10^7),
                                expression(5.0%*%10^7),
                                expression(7.5%*%10^7),
                                expression(10.0%*%10^7)),
                     breaks = c(0,25000000,50000000,75000000,100000000))  +
  theme(panel.grid = element_blank(), 
        panel.background = element_rect(fill = NA), 
        axis.text.y = element_text(color = colors[2],size = 20),
        axis.ticks.y = element_line(color = colors[2]),
        axis.title.x = element_text(size = 23),
        axis.title.y = element_text(color = colors[2],size = 20), 
        axis.line.y = element_line(color = colors[2]), 
        axis.text.x = element_text(size = 20,angle = 45, hjust = 1, vjust = 1)
  )



p3 <- ggplot(dd, aes(Time, Immunity, group = 1)) + 
  geom_line(size=1.5,colour = colors[3]) + 
  scale_y_continuous(limits = c(-6, 20), expand = c(0,0),
                     labels = c(expression(italic(-0.6%*%10^1)),
                                expression(0),
                                expression(0.6%*%10^1),
                                expression(1.2%*%10^1),
                                expression(1.8%*%10^1)),
                     breaks = c(-6,0,6,12,18)) +
  labs(x = "Time(h)", y = expression("Immunity cells(Cells/g)")) +
  theme(panel.grid = element_blank(), 
        panel.background = element_rect(fill = NA), 
        axis.text.y = element_text(color = colors[3],size = 20),
        axis.ticks.y = element_line(color = colors[3]),
        axis.title.x = element_text(size = 23),
        axis.title.y = element_text(color = colors[3],size = 20), 
        axis.line.y = element_line(color = colors[3]), 
        axis.text.x = element_text(size = 20,angle = 45, hjust = 1, vjust = 1)
  )


plot_multi_yaxis(p1, p2, p3)



############################## Immunity existance strong


r = 1
phi = 5*10^-8
eps = 8.2*10^-8
KD = 2.2*10^6
KC = 10^9

beta = 100
phi = 5*10^-8
omg = 1

alpha = 0.97
KI = 2.4*10^7
KN = 10^5
Iomg = 0.05


derivs <- function(t,state, parms) {
  with(as.list(c(state, parms)), {
    
    dB <- r*B*(1-B/KC) - B*phi*P - eps*I*B/(1+B/KD)
    
    dP <- beta*phi*P*B - omg*P
    
    dI <- alpha * I * (1-I/KI) * B/(B + KN) - Iomg*I
    
    list(c(dB,dP,dI))
  })
}

state <- c(B = 10^6,P = 10^7,I = 2.7*10^6)
times <- seq(0, 25, 0.1)

parms <- c(r,phi,eps,KD,KC,
           beta,phi,omg,
           alpha,KI,KN)

yout <- ode(y = state, times = times, func = derivs, parms = parms,method = 'ode45')
#plot(yout)

########################## visualization
library(ggplot2)
options(scipen = 200)

dat = data.frame(Time = rep(yout[,1],3),
                 Value = c(yout[,2],yout[,3],yout[,4]),
                 Type = c(rep('Bacteria',nrow(yout)),
                          rep('Phage',nrow(yout)),
                          rep('Immunity cells',nrow(yout))
                 )
)


dat_b = subset(dat,dat$Type=='Bacteria')
dat_p = subset(dat,dat$Type=='Phage')
dat_i = subset(dat,dat$Type=='Immunity cells')

library(ggplot2)

p4 = ggplot(data=dat_b,aes(x=Time,y=Value))+
  geom_line(size=1.5,color="#FF6600") + theme_classic() + 
  theme(legend.position = 'none') + 
  ylab('CFU/g') + xlab('') +
  scale_y_continuous(labels = c(expression(italic(0)),
                                expression(2.5%*%10^5),
                                expression(5.0%*%10^5),
                                expression(7.5%*%10^5),
                                expression(10.0%*%10^5)),
                     position = "left",
                     expand = c(0,0),
                     breaks = c(0,250000,500000,750000,1000000),
                     limits = c(0,1000000))



p5 = ggplot(data=dat_p,aes(x=Time,y=Value))+
  geom_line(size=1.5,color="#33CC99") + theme_classic() + 
  theme(legend.position = 'none') + 
  ylab('PFU/g') + xlab('') +
  scale_y_continuous(labels = c(expression(italic(0)),
                                expression(2.5%*%10^7),
                                expression(5.0%*%10^7),
                                expression(7.5%*%10^7),
                                expression(10.0%*%10^7)),
                     position = "left",
                     expand = c(0,0),
                     breaks = c(0,25000000,50000000,75000000,100000000),
                     limits = c(0,100000000))


p6 = ggplot(data=dat_i,aes(x=Time,y=Value))+
  geom_line(size=1.5,color="#6633FF") + theme_classic() + 
  theme(legend.position = 'none') + 
  ylab('Cells/g') + xlab('') +
  scale_y_continuous(labels = c(expression(italic(0)),
                                expression(0.5%*%10^7),
                                expression(1.0%*%10^7),
                                expression(1.5%*%10^7),
                                expression(2.0%*%10^7)),
                     position = "left",
                     expand = c(0,0),
                     breaks = c(0,5000000,10000000,15000000,20000000),
                     limits = c(0,20000000))


dd = data.frame(cbind(Time = yout[,1],Bacteria = yout[,2],Phage = yout[,3],Immunity = yout[,4]))

colors <- c('#E83529', '#106EB8', '#F39800')
p1 <- ggplot(dd, aes(Time, Bacteria, group = 1)) + 
  geom_line(size=1.5,colour = colors[1], position = position_nudge(x = -0.2)) + 
  labs(x = "Time(h)", y = "Bacteria(CFU/g)") +
  scale_y_continuous(limits = c(0,1200000), expand = c(0,0),
                     labels = c(expression(italic(0)),
                                expression(2.5%*%10^5),
                                expression(5.0%*%10^5),
                                expression(7.5%*%10^5),
                                expression(10.0%*%10^5)),
                     breaks = c(0,250000,500000,750000,1000000)) +
  theme(panel.grid = element_blank(), 
        panel.background = element_rect(fill = NA), 
        axis.text.y = element_text(color = colors[1],size = 20),
        axis.ticks.y = element_line(color = colors[1]),
        axis.title.x = element_text(size = 23),
        axis.title.y = element_text(color = colors[1],size = 20), 
        axis.line.y = element_line(color = colors[1]),
        axis.line.x = element_line(color = 'black'),
        axis.text.x = element_text(size = 20,angle = 45, hjust = 1, vjust = 1)
  )


p2 <- ggplot(dd, aes(Time, Phage, group = 1)) + 
  geom_line(size=1.5,colour = colors[2]) + 
  labs(x = "Time(h)", y = "Phage(PFU/g)") +
  scale_y_continuous(limits = c(0,120000000), expand = c(0,0),
                     labels = c(expression(italic(0)),
                                expression(2.5%*%10^7),
                                expression(5.0%*%10^7),
                                expression(7.5%*%10^7),
                                expression(10.0%*%10^7)),
                     breaks = c(0,25000000,50000000,75000000,100000000))  +
  theme(panel.grid = element_blank(), 
        panel.background = element_rect(fill = NA), 
        axis.text.y = element_text(color = colors[2],size = 20),
        axis.ticks.y = element_line(color = colors[2]),
        axis.title.x = element_text(size = 23),
        axis.title.y = element_text(color = colors[2],size = 20), 
        axis.line.y = element_line(color = colors[2]),
        axis.line.x = element_line(color = 'black'),
        axis.text.x = element_text(size = 20,angle = 45, hjust = 1, vjust = 1)
  )


p3 <- ggplot(dd, aes(Time, Immunity, group = 1)) + 
  geom_line(size=1,colour = colors[3]) + 
  scale_y_continuous(limits = c(0,22000000), expand = c(0,0),
                     labels = c(expression(italic(0)),
                                expression(0.5%*%10^7),
                                expression(1.0%*%10^7),
                                expression(1.5%*%10^7),
                                expression(2.0%*%10^7)),
                     breaks = c(0,5000000,10000000,15000000,20000000)) +
  labs(x = "Time(h)", y = expression("Immunity cells(Cells/g)")) +
  theme(panel.grid = element_blank(), 
        panel.background = element_rect(fill = NA), 
        axis.text.y = element_text(color = colors[3],size = 20),
        axis.ticks.y = element_line(color = colors[3]),
        axis.title.x = element_text(size = 23),
        axis.title.y = element_text(color = colors[3],size = 20), 
        axis.line.y = element_line(color = colors[3]),
        axis.line.x = element_line(color = 'black'),
        axis.text.x = element_text(size = 20,angle = 45, hjust = 1, vjust = 1)
  )

plot_multi_yaxis(p1, p2, p3)



############################## Immunity existance weak


r = 1
phi = 5*10^-8
eps = 8.2*10^-8
KD = 2.2*10^6
KC = 10^9

beta = 100
phi = 5*10^-8
omg = 1

alpha = 0.05
KI = 2.4*10^7
KN = 10^5
Iomg = 0.005


derivs <- function(t,state, parms) {
  with(as.list(c(state, parms)), {
    
    dB <- r*B*(1-B/KC) - B*phi*P - eps*I*B/(1+B/KD)
    
    dP <- beta*phi*P*B - omg*P
    
    dI <- alpha * I * (1-I/KI) * B/(B + KN) - Iomg*I
    
    list(c(dB,dP,dI))
  })
}

state <- c(B = 10^6,P = 10^7,I = 2.7*10^6)
times <- seq(0, 25, 0.1)

parms <- c(r,phi,eps,KD,KC,
           beta,phi,omg,
           alpha,KI,KN)

yout <- ode(y = state, times = times, func = derivs, parms = parms,method = 'ode45')
#plot(yout)

########################## visualization
library(ggplot2)
options(scipen = 200)

dat = data.frame(Time = rep(yout[,1],3),
                 Value = c(yout[,2],yout[,3],yout[,4]),
                 Type = c(rep('Bacteria',nrow(yout)),
                          rep('Phage',nrow(yout)),
                          rep('Immunity cells',nrow(yout))
                 )
)


dat_b = subset(dat,dat$Type=='Bacteria')
dat_p = subset(dat,dat$Type=='Phage')
dat_i = subset(dat,dat$Type=='Immunity cells')

library(ggplot2)

p7 = ggplot(data=dat_b,aes(x=Time,y=Value))+
  geom_line(size=1.5,color="#FF6600") + theme_classic() + 
  theme(legend.position = 'none') + 
  ylab('CFU/g') + xlab('Time(h)') +
  scale_y_continuous(labels = c(expression(italic(0)),
                                expression(2.5%*%10^5),
                                expression(5.0%*%10^5),
                                expression(7.5%*%10^5),
                                expression(10.0%*%10^5)),
                     position = "left",
                     expand = c(0,0),
                     breaks = c(0,250000,500000,750000,1000000),
                     limits = c(0,1200000))



p8 = ggplot(data=dat_p,aes(x=Time,y=Value))+
  geom_line(size=1.5,color="#33CC99") + theme_classic() + 
  theme(legend.position = 'none') + 
  ylab('PFU/g') + xlab('Time(h)') +
  scale_y_continuous(labels = c(expression(italic(0)),
                                expression(2.5%*%10^7),
                                expression(5.0%*%10^7),
                                expression(7.5%*%10^7),
                                expression(10.0%*%10^7)),
                     position = "left",
                     expand = c(0,0),
                     breaks = c(0,25000000,50000000,75000000,100000000),
                     limits = c(0,110000000))


p9 = ggplot(data=dat_i,aes(x=Time,y=Value))+
  geom_line(size=1.5,color="#6633FF") + theme_classic() + 
  theme(legend.position = 'none') + 
  ylab('Cells/g') + xlab('Time(h)') +
  scale_y_continuous(labels = c(expression(italic(0)),
                                expression(1.0%*%10^6),
                                expression(2.0%*%10^6),
                                expression(3.0%*%10^6),
                                expression(4.0%*%10^6)),
                     position = "left",
                     expand = c(0,0),
                     breaks = c(0,1000000,2000000,3000000,4000000),
                     limits = c(0,4000000))


dd = data.frame(cbind(Time = yout[,1],Bacteria = yout[,2],Phage = yout[,3],Immunity = yout[,4]))

colors <- c('#E83529', '#106EB8', '#F39800')
p1 <- ggplot(dd, aes(Time, Bacteria, group = 1)) + 
  geom_line(size=1.5,colour = colors[1], position = position_nudge(x = -0.2)) + 
  labs(x = "Time(h)", y = "Bacteria(CFU/g)") +
  scale_y_continuous(limits = c(0,1200000), expand = c(0,0),
                     labels = c(expression(italic(0)),
                                expression(2.5%*%10^5),
                                expression(5.0%*%10^5),
                                expression(7.5%*%10^5),
                                expression(10.0%*%10^5)),
                     breaks = c(0,250000,500000,750000,1000000)) +
  theme(panel.grid = element_blank(), 
        panel.background = element_rect(fill = NA), 
        axis.text.y = element_text(color = colors[1],size = 20),
        axis.ticks.y = element_line(color = colors[1]),
        axis.title.x = element_text(size = 23),
        axis.title.y = element_text(color = colors[1],size = 20), 
        axis.line.y = element_line(color = colors[1]),
        axis.line.x = element_line(color = 'black'),
        axis.text.x = element_text(size = 20,angle = 45, hjust = 1, vjust = 1)
  )



p2 <- ggplot(dd, aes(Time, Phage, group = 1)) + 
  geom_line(size=1.5,colour = colors[2]) + 
  labs(x = "Time(h)", y = "Phage(PFU/g)") +
  scale_y_continuous(limits = c(0,120000000), expand = c(0,0),
                     labels = c(expression(italic(0)),
                                expression(2.5%*%10^7),
                                expression(5.0%*%10^7),
                                expression(7.5%*%10^7),
                                expression(10.0%*%10^7)),
                     breaks = c(0,25000000,50000000,75000000,100000000))  +
  theme(panel.grid = element_blank(), 
        panel.background = element_rect(fill = NA), 
        axis.text.y = element_text(color = colors[2],size = 20),
        axis.ticks.y = element_line(color = colors[2]),
        axis.title.x = element_text(size = 23),
        axis.title.y = element_text(color = colors[2],size = 20), 
        axis.line.y = element_line(color = colors[2]),
        axis.text.x = element_text(size = 20,angle = 45, hjust = 1, vjust = 1)
  )


p3 <- ggplot(dd, aes(Time, Immunity, group = 1)) + 
  geom_line(size=1.5,colour = colors[3]) + 
  scale_y_continuous(limits = c(0,22000000), expand = c(0,0),
                     labels = c(expression(italic(0)),
                                expression(0.5%*%10^7),
                                expression(1.0%*%10^7),
                                expression(1.5%*%10^7),
                                expression(2.0%*%10^7)),
                     breaks = c(0,5000000,10000000,15000000,20000000)) +
  labs(x = "Time(h)", y = expression("Immunity cells(Cells/g)")) +
  theme(panel.grid = element_blank(), 
        panel.background = element_rect(fill = NA), 
        axis.text.y = element_text(color = colors[3],size = 20),
        axis.ticks.y = element_line(color = colors[3]),
        axis.title.x = element_text(size = 23),
        axis.title.y = element_text(color = colors[3],size = 20), 
        axis.line.y = element_line(color = colors[3]),
        axis.text.x = element_text(size = 20,angle = 45, hjust = 1, vjust = 1)
  )



plot_multi_yaxis(p1, p2, p3)



################################
############################## Immunity existance weak without phage


r = 1
phi = 0
eps = 5*10^-8
KD = 2.2*10^6
KC = 10^9

beta = 100
phi = 0
omg = 0

alpha = 0.97
KI = 2.4*10^7
KN = 10^5
Iomg = 0.5


derivs <- function(t,state, parms) {
  with(as.list(c(state, parms)), {
    
    dB <- r*B*(1-B/KC) - B*phi*P - eps*I*B/(1+B/KD)
    
    dP <- beta*phi*P*B - omg*P
    
    dI <- alpha * I * (1-I/KI) * B/(B + KN) - Iomg*I
    
    list(c(dB,dP,dI))
  })
}

state <- c(B = 10^6,P = 0,I = 2.7*10^6)
times <- seq(0, 25, 0.1)

parms <- c(r,phi,eps,KD,KC,
           beta,phi,omg,
           alpha,KI,KN)

yout <- ode(y = state, times = times, func = derivs, parms = parms,method = 'ode45')
plot(yout)


########################## visualization
library(ggplot2)
options(scipen = 200)

dat = data.frame(Time = rep(yout[,1],3),
                 Value = c(yout[,2],yout[,3],yout[,4]),
                 Type = c(rep('Bacteria',nrow(yout)),
                          rep('Phage',nrow(yout)),
                          rep('Immunity cells',nrow(yout))
                 )
)


dat_b = subset(dat,dat$Type=='Bacteria')
dat_p = subset(dat,dat$Type=='Phage')
dat_i = subset(dat,dat$Type=='Immunity cells')


library(ggplot2)

p10 = ggplot(data=dat_b,aes(x=Time,y=Value))+
  geom_line(size=1.5,color="#FF6600") + theme_classic() + 
  theme(legend.position = 'none') + 
  ylab('CFU/g') + xlab('Time(h)') +
  scale_y_continuous(labels = c(expression(italic(0)),
                                expression(2.5%*%10^8),
                                expression(5.0%*%10^8),
                                expression(7.5%*%10^8),
                                expression(10.0%*%10^8)),
                     position = "left",
                     expand = c(0,0),
                     breaks = c(0,250000000,500000000,750000000,1000000000),
                     limits = c(0,1100000000))



p11 = ggplot(data=dat_p,aes(x=Time,y=Value))+
  geom_line(size=1.5,color="#33CC99") + theme_classic() + 
  theme(legend.position = 'none') + 
  ylab('PFU/g') + xlab('Time(h)')



p12 = ggplot(data=dat_i,aes(x=Time,y=Value))+
  geom_line(size=1.5,color="#6633FF") + theme_classic() + 
  theme(legend.position = 'none') + 
  ylab('Cells/g') + xlab('Time(h)')
  scale_y_continuous(labels = c(expression(italic(0)),
                                expression(2.5%*%10^6),
                                expression(5.0%*%10^6),
                                expression(7.5%*%10^6),
                                expression(10.0%*%10^6)),
                     position = "left",
                     expand = c(0,0),
                     breaks = c(0,2500000,5000000,7500000,10000000),
                     limits = c(0,12500000))


dd = data.frame(cbind(Time = yout[,1],Bacteria = yout[,2],Phage = yout[,3],Immunity = yout[,4]))

colors <- c('#E83529', '#106EB8', '#F39800')
p1 <- ggplot(dd, aes(Time, Bacteria, group = 1)) + 
  geom_line(size=1,colour = colors[1], position = position_nudge(x = -0.2)) + 
  labs(x = "Time(h)", y = "Bacteria(CFU/g)") +
  scale_y_continuous(limits = c(0,1100000000), expand = c(0,0),
                     llabels = c(expression(italic(0)),
                                 expression(2.5%*%10^8),
                                 expression(5.0%*%10^8),
                                 expression(7.5%*%10^8),
                                 expression(10.0%*%10^8)),
                     breaks = c(0,250000000,500000000,750000000,1000000000)) +
  theme(panel.grid = element_blank(), 
        panel.background = element_rect(fill = NA), 
        axis.text.y = element_text(color = colors[1],size = 15), 
        axis.ticks.y = element_line(color = colors[1]), 
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(color = colors[1],size = 20), 
        axis.line.y = element_line(color = colors[1]), 
        axis.line.x = element_line(color = 'black'),
        axis.text.x = element_text(size = 15,angle = 45, hjust = 1, vjust = 1)
  )


p2 <- ggplot(dd, aes(Time, Phage, group = 1)) + 
  geom_line(size=1,colour = colors[2]) + 
  labs(x = "Time(h)", y = "Phage(PFU/g)") +
  scale_y_continuous(limits = c(-5,20), expand = c(0,0))  +
  theme(panel.grid = element_blank(), 
        panel.background = element_rect(fill = NA), 
        axis.text.y = element_text(color = colors[2],size = 15), 
        axis.ticks.y = element_line(color = colors[2]),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(color = colors[2],size = 20), 
        axis.line.y = element_line(color = colors[2]), 
        axis.text.x = element_text(size = 15,angle = 45, hjust = 1, vjust = 1)
  )

p3 <- ggplot(dd, aes(Time, Immunity, group = 1)) + 
  geom_line(size=1,colour = colors[3]) + 
  scale_y_continuous(limits = c(0,12500000), expand = c(0,0),
                     labels = c(expression(italic(0)),
                                expression(2.5%*%10^6),
                                expression(5.0%*%10^6),
                                expression(7.5%*%10^6),
                                expression(10.0%*%10^6)),
                     breaks = c(0,2500000,5000000,7500000,10000000)) +
  labs(x = "Time(h)", y = expression("Immunity cells(Cells/g)")) +
  theme(panel.grid = element_blank(), 
        panel.background = element_rect(fill = NA), 
        axis.text.y = element_text(color = colors[3],size = 15), 
        axis.ticks.y = element_line(color = colors[3]),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(color = colors[3],size = 20), 
        axis.line.y = element_line(color = colors[3]), 
        axis.text.x = element_text(size = 15,angle = 45, hjust = 1, vjust = 1)
  )


plot_multi_yaxis(p1, p2, p3)


  
  
############################## Immunity existance strong without phage


r = 1
phi = 0
eps = 16.4*10^-8
KD = 2.2*10^6
KC = 10^9

beta = 100
phi = 0
omg = 0

alpha = 0.97
KI = 2.4*10^7
KN = 10^5
Iomg = 0.05


derivs <- function(t,state, parms) {
  with(as.list(c(state, parms)), {
    
    dB <- r*B*(1-B/KC) - B*phi*P - eps*I*B/(1+B/KD)
    
    dP <- beta*phi*P*B - omg*P
    
    dI <- alpha * I * (1-I/KI) * B/(B + KN) - Iomg*I
    
    list(c(dB,dP,dI))
  })
}

state <- c(B = 10^6,P = 0,I = 2.7*10^6)
times <- seq(0, 25, 0.1)

parms <- c(r,phi,eps,KD,KC,
           beta,phi,omg,
           alpha,KI,KN)

yout <- ode(y = state, times = times, func = derivs, parms = parms,method = 'ode45')
plot(yout)

########################## visualization
library(ggplot2)
options(scipen = 200)

dat = data.frame(Time = rep(yout[,1],3),
                 Value = c(yout[,2],yout[,3],yout[,4]),
                 Type = c(rep('Bacteria',nrow(yout)),
                          rep('Phage',nrow(yout)),
                          rep('Immunity cells',nrow(yout))
                 )
)


dat_b = subset(dat,dat$Type=='Bacteria')
dat_p = subset(dat,dat$Type=='Phage')
dat_i = subset(dat,dat$Type=='Immunity cells')


library(ggplot2)

p13 = ggplot(data=dat_b,aes(x=Time,y=Value))+
  geom_line(size=1.5,color="#FF6600") + theme_classic() + 
  theme(legend.position = 'none') + 
  ylab('CFU/g') + xlab('Time(h)') +
  scale_y_continuous(labels = c(expression(italic(0)),
                                expression(1.0%*%10^6),
                                expression(2.0%*%10^6),
                                expression(3.0%*%10^6)),
                     position = "left",
                     expand = c(0,0),
                     breaks = c(0,1000000,2000000,3000000),
                     limits = c(0,3300000))



p14 = ggplot(data=dat_p,aes(x=Time,y=Value))+
  geom_line(size=1.5,color="#33CC99") + theme_classic() + 
  theme(legend.position = 'none') + 
  ylab('PFU/g') + xlab('Time(h)')



p15 = ggplot(data=dat_i,aes(x=Time,y=Value))+
  geom_line(size=1.5,color="#6633FF") + theme_classic() + 
  theme(legend.position = 'none') + 
  ylab('Cells/g') + xlab('Time(h)') + 
  scale_y_continuous(labels = c(expression(italic(0)),
                                expression(5.0%*%10^6),
                                expression(10.0%*%10^6),
                                expression(15.0%*%10^6),
                                expression(20.0%*%10^6)),
                     position = "left",
                     expand = c(0,0),
                     breaks = c(0,5000000,10000000,15000000,20000000),
                     limits = c(0,23000000))



dd = data.frame(cbind(Time = yout[,1],Bacteria = yout[,2],Phage = yout[,3],Immunity = yout[,4]))

colors <- c('#E83529', '#106EB8', '#F39800')
p1 <- ggplot(dd, aes(Time, Bacteria, group = 1)) + 
  geom_line(size=1,colour = colors[1], position = position_nudge(x = -0.2)) + 
  labs(x = "Time(h)", y = "Bacteria(CFU/g)") +
  scale_y_continuous(limits = c(0,3300000), expand = c(0,0),
                     labels = c(expression(italic(0)),
                                expression(1.0%*%10^6),
                                expression(2.0%*%10^6),
                                expression(3.0%*%10^6)),
                     breaks = c(0,1000000,2000000,3000000)) +
  theme(panel.grid = element_blank(), 
        panel.background = element_rect(fill = NA), 
        axis.text.y = element_text(color = colors[1],size = 15), 
        axis.ticks.y = element_line(color = colors[1]), 
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(color = colors[1],size = 20), 
        axis.line.y = element_line(color = colors[1]), 
        axis.line.x = element_line(color = 'black'),
        axis.text.x = element_text(size = 15,angle = 45, hjust = 1, vjust = 1)
  )


p2 <- ggplot(dd, aes(Time, Phage, group = 1)) + 
  geom_line(size=1,colour = colors[2]) + 
  labs(x = "Time(h)", y = "Phage(PFU/g)") +
  scale_y_continuous(limits = c(-5,20), expand = c(0,0))  +
  theme(panel.grid = element_blank(), 
        panel.background = element_rect(fill = NA), 
        axis.text.y = element_text(color = colors[2],size = 15), 
        axis.ticks.y = element_line(color = colors[2]),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(color = colors[2],size = 20), 
        axis.line.y = element_line(color = colors[2]), 
        axis.text.x = element_text(size = 15,angle = 45, hjust = 1, vjust = 1)
  )


p3 <- ggplot(dd, aes(Time, Immunity, group = 1)) + 
  geom_line(size=1,colour = colors[3]) + 
  scale_y_continuous(limits = c(0,23000000), expand = c(0,0),
                     labels = c(expression(italic(0)),
                                expression(1.0%*%10^6),
                                expression(2.0%*%10^6),
                                expression(3.0%*%10^6)),
                     breaks = c(0,5000000,10000000,15000000,20000000)) +
  labs(x = "Time(h)", y = expression("Immunity cells(Cells/g)")) +
  theme(panel.grid = element_blank(), 
        panel.background = element_rect(fill = NA), 
        axis.text.y = element_text(color = colors[3],size = 15), 
        axis.ticks.y = element_line(color = colors[3]),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(color = colors[3],size = 20), 
        axis.line.y = element_line(color = colors[3]), 
        axis.text.x = element_text(size = 15,angle = 45, hjust = 1, vjust = 1)
  )

plot_multi_yaxis(p1, p2, p3)


##############################
library(patchwork)

(p1 + p2 + p3) / (p4 + p5 + p6) / (p7 + p8 + p9) / (p10 + p11 + p12) / (p13 + p14 + p15)


(p1 + p2 + p3) / (p4 + p5 + p6) / (p7 + p8 + p9)

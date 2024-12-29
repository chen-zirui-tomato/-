# 设置标题和轴标签
set title "三次B样条插值"
set xlabel "x"
set ylabel "y"

# 设置输出文件格式（例如 PNG）并输出第一张图
set terminal pngcairo size 800,600
set output 'figure/splineA.png'  # 第一张图

# 绘制数据
plot 'bin/testA.dat' with linespoints title '插值曲线'
# 关闭输出
unset output

# 设置输出文件格式并输出第二张图
set output 'figure/splineB.png'  # 第二张图

# 绘制数据
plot 'bin/testB.dat' with linespoints title '插值曲线'
# 关闭输出
unset output

#################################################################

# 设置标题和轴标签
set title "心型线B样条拟合"
set xlabel "x"
set ylabel "y"


# 设置输出文件格式（例如 PNG）并输出图
set terminal pngcairo size 800,600
set output 'figure/splineE1.png'  # 合并后的图

# 设置颜色
set style line 1 linecolor rgb 'blue' linetype 1 linewidth 2 # 使用蓝色作为线条颜色

# 绘制数据，使用同一种颜色，只显示一个标题
plot 'bin/testE1.dat' with linespoints linestyle 1 title 'E1 插值曲线', \
     'bin/testE2.dat' with linespoints linestyle 1 title ''

# 关闭输出
unset output

#############################################################

# 设置输出文件格式（例如 PNG）并输出图
set terminal pngcairo size 800,600
set output 'figure/splineE2.png'  # 合并后的图

# 设置颜色
set style line 1 linecolor rgb 'blue' linetype 1 linewidth 2 # 使用蓝色作为线条颜色

# 绘制数据，使用同一种颜色，只显示一个标题
plot 'bin/testE3.dat' with linespoints linestyle 1 title 'E1 插值曲线', \
     'bin/testE4.dat' with linespoints linestyle 1 title ''

# 关闭输出
unset output

#############################################################

# 设置标题和轴标签
set title "螺线B样条拟合"
set xlabel "x"
set ylabel "y"

# 设置输出文件格式（例如 PNG）并输出图
set terminal pngcairo size 800,600
set output 'figure/splineE3.png'  

# 设置颜色
set style line 1 linecolor rgb 'blue' linetype 1 linewidth 2 # 使用蓝色作为线条颜色

# 绘制数据
plot 'bin/testE5.dat' with linespoints title '插值曲线'

# 关闭输出
unset output

###############################################################
# 设置输出文件格式（例如 PNG）并输出图
set terminal pngcairo size 800,600
set output 'figure/splineE4.png'  

# 设置颜色
set style line 1 linecolor rgb 'blue' linetype 1 linewidth 2 # 使用蓝色作为线条颜色

# 绘制数据
plot 'bin/testE6.dat' with linespoints title '插值曲线'

# 关闭输出
unset output


##################################################################
set terminal pngcairo size 800,600
set output 'figure/threeDplot.png'

set style data points
set pointsize 1.5
set title "三维数据点"
set xlabel "X轴"
set ylabel "Y轴"
set zlabel "Z轴"

splot 'bin/testE7.dat' with points title '数据点'

unset output

################################################################

set output 'figure/splineM.png'  # 第二张图

# 绘制数据
plot 'bin/testM.dat' with linespoints title '插值曲线'
# 关闭输出
unset output

################################################################

set output 'figure/splineAR.png'  # 第二张图

# 绘制数据
plot 'bin/testAR.dat' with linespoints title '插值曲线'
# 关闭输出
unset output
CXX = g++

# 编译选项
CXXFLAGS = -std=c++11 -Wall -Wextra

# 源文件
# SRCS = ProblemF.cpp
SRCS = $(wildcard Problem*.cpp)

# 目标文件（提取文件名去掉前缀 "Problem"）
TARGETS = $(SRCS:Problem%.cpp=%)

# LaTeX 文件
TEXFILE = report.tex
PDF = $(TEXFILE:.tex=.pdf)

# 默认规则
all: $(TARGETS) $(PDF)

# 编译目标
%: Problem%.cpp
	$(CXX) $(CXXFLAGS) -o $@ $<

# 编译 LaTeX 文件
$(PDF): $(TEXFILE)
	xelatex $(TEXFILE)

# 清理目标
clean:
	rm -f $(TARGETS) $(PDF) *.aux *.log *.fls *.fdb_latexmk *.synctex* output*.dat

.PHONY: all clean

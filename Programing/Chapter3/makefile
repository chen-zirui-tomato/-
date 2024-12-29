# 定义编译器
CXX = g++

# 定义编译选项
CXXFLAGS = -std=c++11 -Wall -Wextra -I /usr/include/eigen3

# 添加 Boost 库的链接选项
LDFLAGS = -lboost_system -lboost_filesystem 

# 定义源文件路径
SRC_DIR = src
INC_DIR = include
TAR_DIR = bin
PIC_DIR = figure

# 定义源文件
SRCS = $(wildcard $(SRC_DIR)/Problem*.cpp)

# 定义目标文件
OBJS = $(SRCS:.cpp=.o)

# 定义目标程序名称，并指定生成路径
TARGETS = $(patsubst $(SRC_DIR)/Problem%.cpp, $(TAR_DIR)/Problem% , $(SRCS))

# 默认目标
all: run

# 链接目标
$(TAR_DIR)/Problem%: $(SRC_DIR)/Problem%.o
	$(CXX) $(CXXFLAGS) -o $@ $<

# 生成目标文件
$(SRC_DIR)/%.o: $(SRC_DIR)/%.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# LaTeX 文件路径
DOC_DIR = doc
TEXFILE = $(DOC_DIR)/design.tex
PDF = $(TEXFILE:.tex=.pdf)

# 默认规则
all: $(TARGETS) $(PDF)

# 编译 LaTeX 文件
$(PDF): $(TEXFILE)
	cd $(DOC_DIR) && xelatex $(notdir $(TEXFILE))


# 运行规则
run: $(TARGETS)
	@for target in $(TARGETS); do \
		./$$target; \
	done
	gnuplot $(SRC_DIR)/plot.gp

# 清理生成的文件
clean:
	rm -f $(OBJS) $(TAR_DIR)/Problem* $(SRC_DIR)/test.dat $(SRC_DIR)/*.png 
	rm -f $(DOC_DIR)/*.aux $(DOC_DIR)/*.log $(DOC_DIR)/*.pdf $(DOC_DIR)/*.out

# 伪目标
.PHONY: all run clean

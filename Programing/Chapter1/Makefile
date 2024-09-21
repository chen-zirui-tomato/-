# Makefile for compiling B.cpp

# 编译器
CXX = g++

# 编译选项
CXXFLAGS = -std=c++11 -Wall -Wextra

# 目标文件
TARGET = B

# 源文件
SRCS = ProblemF.cpp

# 头文件
HEADERS = Function.hpp EquationSolver.hpp

# 默认规则
all: $(TARGET)

# 编译目标
$(TARGET): $(SRCS) $(HEADERS)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(SRCS)

# 清理目标
clean:
	rm -f $(TARGET)

.PHONY: all clean

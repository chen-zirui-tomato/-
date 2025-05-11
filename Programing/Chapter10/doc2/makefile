# LaTeX 文件
TEXFILE = report.tex
PDF = $(TEXFILE:.tex=.pdf)

# 默认规则
all: $(PDF)

# 编译 LaTeX 文件
$(PDF): $(TEXFILE)
	xelatex $(TEXFILE)
	xelatex $(TEXFILE)

# 清理生成的临时文件
clean:
	rm -f $(PDF) $(TEXFILE:.tex=.aux) $(TEXFILE:.tex=.log) $(TEXFILE:.tex=.out)

# 伪目标
.PHONY: all clean

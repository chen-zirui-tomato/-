# LaTeX 文件
TEXFILE = report.tex
PDF = $(TEXFILE:.tex=.pdf)

# 默认规则
all: $(TARGETS) $(PDF)

# 编译 LaTeX 文件
$(PDF): $(TEXFILE)
	xelatex $(TEXFILE)

# 清理目标
clean:
	rm -f $(TARGETS) $(PDF) *.aux *.log *.fls *.fdb_latexmk *.synctex* output*.dat

.PHONY: all clean

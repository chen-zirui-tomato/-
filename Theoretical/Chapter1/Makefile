# 定义变量
TEXFILE = Homework.tex
PDFFILE = $(TEXFILE:.tex=.pdf)
BIBFILE = citation.bib

# 默认目标
all: $(PDFFILE)

# 生成 PDF
$(PDFFILE): $(TEXFILE) $(BIBFILE)
	pdflatex $(TEXFILE)
	biber $(TEXFILE:.tex=)
	pdflatex $(TEXFILE)
	pdflatex $(TEXFILE)

# 清理生成的文件
clean:
	rm -f *.aux *.log *.bbl *.bcf *.blg *.out *.run.xml *.pdf

.PHONY: all clean

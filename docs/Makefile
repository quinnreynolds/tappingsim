SPHINXBUILD = python3 -m sphinx
SOURCEDIR   = .
BUILDDIR    = _build

help:
	@$(SPHINXBUILD) -M help $(SOURCEDIR) $(BUILDDIR) $(SPHINXOPTS) $(O)

.PHONY: help Makefile html clean

html: notebooks
	@$(SPHINXBUILD) -M html $(SOURCEDIR) $(BUILDDIR) $(SPHINXOPTS) $(O)

notebooks:
	@cp -r ../notebooks ./notebooks

clean:
	@rm -rf $(BUILDDIR) notebooks

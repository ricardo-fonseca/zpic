BUILD_LIB = ../lib
DOC_HTML = ../../documentation/python/api
DOC_MODS = em1d em2d em1ds em2ds es1d

all:
	python3 setup.py build_ext --build-lib=$(BUILD_LIB)

docs : all
	@echo "(*info*) PEP-224 Warnings can be safely ignored"
	pdoc3 --html --force --output-dir $(DOC_HTML) $(DOC_MODS)

clean:
	rm -rf build *.c $(BUILD_LIB)/*.so

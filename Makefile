# This file is part of schema-salad,
# https://github.com/common-workflow-language/schema-salad/, and is
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#
# Contact: common-workflow-language@googlegroups.com

# make format to fix most python formatting errors
# make pylint to check Python code for enhanced compliance including naming
#  and documentation
# make coverage-report to check coverage of the python scripts by the tests

MODULE=wdl2cwl
PACKAGE=wdl2cwl

# `SHELL=bash` doesn't work for some, so don't use BASH-isms like
# `[[` conditional expressions.
PYSOURCES=$(wildcard ${MODULE}/**.py ${MODULE}/avro/*.py ${MODULE}/tests/*.py) setup.py
DEVPKGS=diff_cover pylint pep257 pydocstyle flake8 tox tox-pyenv \
	isort wheel autoflake pyupgrade bandit -rlint-requirements.txt \
	-rtest-requirements.txt -rmypy_requirements.txt
COVBASE=coverage run --append
PYTEST_EXTRA ?= -rs

# Updating the Major & Minor version below?
# Don't forget to update setup.py as well
VERSION=7.1.$(shell date +%Y%m%d%H%M%S --utc --date=`git log --first-parent \
	--max-count=1 --format=format:%cI`)

## all         : default task
all: dev

## help        : print this help message and exit
help: Makefile
	@sed -n 's/^##//p' $<

## cleanup                : shortcut for "make sort_imports format flake8 diff_pydocstyle_report"
cleanup: sort_imports format flake8 diff_pydocstyle_report

## install-dep : install most of the development dependencies via pip
install-dep: install-dependencies

install-dependencies: FORCE
	pip install --upgrade $(DEVPKGS)
	pip install -r requirements.txt

## install     : install the ${MODULE} module and schema-salad-tool
install: FORCE
	pip install .

## dev     : install the ${MODULE} module in dev mode
dev: install-dep
	pip install -e .

## dist        : create a module package for distribution
dist: dist/${MODULE}-$(VERSION).tar.gz

dist/${MODULE}-$(VERSION).tar.gz: $(SOURCES)
	python setup.py sdist bdist_wheel

## docs	       : make the docs
docs: FORCE
	cd docs && $(MAKE) html

## clean       : clean up all temporary / machine-generated files
clean: FORCE
	rm -f ${MODILE}/*.pyc tests/*.pyc
	python setup.py clean --all || true
	rm -Rf .coverage coverage.xml
	rm -f diff-cover.html

# Linting and code style related targets
## sorting imports using isort: https://github.com/timothycrosley/isort
sort_imports: $(filter-out wdl2cwl/WdlV1_1%,$(PYSOURCES))
	isort $^

remove_unused_imports: $(filter-out wdl2cwl/WdlV1_1%,$(PYSOURCES))
	autoflake --in-place --remove-all-unused-imports $^

pep257: pydocstyle
## pydocstyle      : check Python code style
pydocstyle: $(PYSOURCES)
	pydocstyle --add-ignore=D100,D101,D102,D103 $^ || true

pydocstyle_report.txt: $(PYSOURCES)
	pydocstyle setup.py $^ > $@ 2>&1 || true

diff_pydocstyle_report: pydocstyle_report.txt
	diff-quality --compare-branch=main --violations=pydocstyle --fail-under=100 $^

## format      : check/fix all code indentation and formatting (runs black)
format: $(filter-out wdl2cwl/WdlV1_1%,$(PYSOURCES))
	black $^

format-check: $(filter-out wdl2cwl/WdlV1_1%,$(PYSOURCES))
	black --diff --check $^

## pylint      : run static code analysis on Python code
pylint: $(PYSOURCES)
	pylint --msg-template="{path}:{line}: [{msg_id}({symbol}), {obj}] {msg}" \
                $^ -j0|| true

pylint_report.txt: $(PYSOURCES)
	pylint --msg-template="{path}:{line}: [{msg_id}({symbol}), {obj}] {msg}" \
		$^ -j0> $@ || true

diff_pylint_report: pylint_report.txt
	diff-quality --violations=pylint pylint_report.txt

.coverage:
	pytest --cov --cov-config=.coveragerc --cov-report= -n auto ${PYTEST_EXTRA}

coverage.xml: .coverage
	coverage xml

coverage.html: htmlcov/index.html

htmlcov/index.html: .coverage
	coverage html
	@echo Test coverage of the Python code is now in htmlcov/index.html

coverage-report: .coverage
	coverage report

diff-cover: coverage.xml
	diff-cover $^

diff-cover.html: coverage.xml
	diff-cover $^ --html-report $@

## test        : run the ${MODULE} test suite
test: $(PYSOURCES)
	python -m pytest ${PYTEST_EXTRA}

vpath %.wdl wdl2cwl/tests/wdl_files

wdl2cwl/tests/cwl_files/%.cwl: %.wdl wdl2cwl/main.py
	wdl2cwl $< --output $@

validate-cwl:
	$(foreach file,$(wildcard wdl2cwl/tests/cwl_files/*.cwl),cwltool --validate $(file) &&) true

# See https://github.com/common-workflow-lab/wdl2cwl_test_cache
cwltest:
	cd wdl2cwl/tests && cwltest --timeout 1200 --verbose --junit-verbose --tool cwltool --test cwl_tests.yaml -- --cachedir cache

cwltest_bigdata:
	cd wdl2cwl/tests && cwltest --timeout 1200 --verbose --junit-verbose --tool cwltool --test cwl_tests_bigdata.yaml -- --cachedir cache 

cwltest_bigdata_bigmem:
	cd wdl2cwl/tests && cwltest --timeout 1200 --verbose --junit-verbose --tool cwltool --test cwl_tests_bigdata_bigmem.yaml -- --cachedir cache

## testcov     : run the ${MODULE} test suite and collect coverage
testcov: $(PYSOURCES)
	pytest --cov ${PYTEST_EXTRA}

sloccount.sc: $(PYSOURCES) Makefile
	sloccount --duplicates --wide --details $^ > $@

## sloccount   : count lines of code
sloccount: $(PYSOURCES) Makefile
	sloccount $^

list-author-emails:
	@echo 'name, E-Mail Address'
	@git log --format='%aN,%aE' | sort -u | grep -v 'root'

mypy3: mypy
mypy: $(filter-out setup.py wdl2cwl/WdlV1_1%,$(PYSOURCES))
	mypy --exclude '\.tox.*' $^

pyupgrade: $(filter-out wdl2cwl/WdlV1_1%,$(PYSOURCES))
	pyupgrade --exit-zero-even-if-changed --py37-plus $^

release-test: FORCE
	git diff-index --quiet HEAD -- || ( echo You have uncommitted changes, please commit them and try again; false )
	./release-test.sh

release: release-test
	. testenv2/bin/activate && \
		python testenv2/src/${PACKAGE}/setup.py sdist bdist_wheel
	. testenv2/bin/activate && \
		pip install twine && \
		twine upload testenv2/src/${PACKAGE}/dist/* && \
		git tag ${VERSION} && git push --tags

flake8: $(filter-out wdl2cwl/WdlV1_1%,$(PYSOURCES))
	flake8 $^

bandit: $(filter-out wdl2cwl/tests% wdl2cwl/WdlV1_1%,$(PYSOURCES))
	bandit $^

FORCE:

# Use this to print the value of a Makefile variable
# Example `make print-VERSION`
# From https://www.cmcrossroads.com/article/printing-value-makefile-variable
print-%  : ; @echo $* = $($*)

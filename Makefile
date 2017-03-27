
PYTHON ?= python2.7

# For these tests, you need to check out testdata from a different
# repository:
#
# git clone https://github.com/legacysurvey/obsbot-testdata.git testdata
#

test:
	$(PYTHON) obsdb/manage.py test -p test_copilot.py
	$(PYTHON) obsdb/manage.py test -p test_mosbot.py

test_mosbot:
	$(PYTHON) obsdb/manage.py test -p test_mosbot.py

test_copilot:
	$(PYTHON) obsdb/manage.py test -p test_copilot.py

test_copilot_one:
	$(PYTHON) obsdb/manage.py test test_copilot.TestCopilot.test_new_image

test_focus:
	$(PYTHON) obsdb/manage.py test test_copilot.TestCopilot.test_focus

test_decbot:
	$(PYTHON) test_decbot.py
	$(PYTHON) obsdb/manage.py test -p test_decbot_2.py

coverage_decbot:
	coverage erase
	coverage run test_decbot.py
	coverage run -a obsdb/manage.py test -p test_decbot_2.py
	coverage report
	coverage html
	@echo
	@echo "Now you might want to:"
	@echo "  open coverage_html_report/index.html"
	@echo

.PHONY: coverage_decbot

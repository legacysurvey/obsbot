
PYTHON ?= python2.7

test:
	$(PYTHON) obsdb/manage.py test -p test_copilot.py
	$(PYTHON) test_mosbot.py

test_copilot:
	$(PYTHON) obsdb/manage.py test -p test_copilot.py

test_copilot_one:
	$(PYTHON) obsdb/manage.py test test_copilot.TestCopilot.test_new_image

test_focus:
	$(PYTHON) obsdb/manage.py test test_copilot.TestCopilot.test_focus

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

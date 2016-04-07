
PYTHON ?= python2.7

test:
	$(PYTHON) obsdb/manage.py test -p test_copilot.py
	$(PYTHON) test_mosbot.py

test_copilot:
	$(PYTHON) obsdb/manage.py test -p test_copilot.py

test_focus:
	$(PYTHON) obsdb/manage.py test test_copilot.TestCopilot.test_focus

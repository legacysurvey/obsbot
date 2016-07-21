
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

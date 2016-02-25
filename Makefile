
test:
	python obsdb/manage.py test -p test_copilot.py
	python test_mosbot.py


test_focus:
	python obsdb/manage.py test test_copilot.TestCopilot.test_focus

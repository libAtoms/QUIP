import unittest, os.path, glob


test_files = glob.glob(os.path.join(os.path.dirname(__file__), 'test*.py'))
test_mods = [ os.path.splitext(os.path.basename(t))[0] for t in test_files ]

suite = unittest.TestLoader().loadTestsFromNames(test_mods)
unittest.TextTestRunner().run(suite)
   

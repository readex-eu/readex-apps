import unittest
import textwrap

class TestToolOutput(unittest.TestCase):
    
    def setUp(self):
        with open("tool_output") as f:
            self.output = f.read()

    def test_numthreads_2(self):
        expected_output = textwrap.dedent("""\
            [OpenMPTP][NUMTHREADS]: Get old setting
            [OpenMPTP][NUMTHREADS]: Old Setting = 24
            [OpenMPTP][NUMTHREADS]: Set new setting
            [OpenMPTP][NUMTHREADS]: New Setting = 2""")

        self.assertIn(expected_output, self.output)

    def test_numthreads_1(self):
        expected_output = textwrap.dedent("""\
            [OpenMPTP][NUMTHREADS]: Get old setting
            [OpenMPTP][NUMTHREADS]: Old Setting = 2
            [OpenMPTP][NUMTHREADS]: Set new setting
            [OpenMPTP][NUMTHREADS]: New Setting = 1""")

        self.assertIn(expected_output, self.output)
        
    def test_no_exceptions(self):
        not_expected_output = textwrap.dedent("""\
            Uncaught exception with message:""")
        self.assertNotIn(not_expected_output, self.output)


if __name__ == '__main__':
    unittest.main()

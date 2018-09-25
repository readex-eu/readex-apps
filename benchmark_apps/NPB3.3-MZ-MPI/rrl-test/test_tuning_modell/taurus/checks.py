import unittest
import textwrap

class TestToolOutput(unittest.TestCase):
    
    def setUp(self):
        with open("tool_output") as f:
            self.output = f.read()

    def test_omp_shedule_type_2_to_1(self):
        expected_output = textwrap.dedent("""\
            [OpenMPTP][SCHEDULE_TYPE]: setting scheduling type
            [OpenMPTP][SCHEDULE_TYPE]: Old kind = 2 Old chunk size = 1
            [OpenMPTP][SCHEDULE_TYPE]: New kind = 1 New chunk size = 1""")
        self.assertIn(expected_output, self.output)

    def test_omp_chunk_size_1_to_4(self):
        expected_output = textwrap.dedent("""\
            [OpenMPTP][SCHEDULE_CHUNK_SIZE]: setting scheduling chunk size
            [OpenMPTP][SCHEDULE_CHUNK_SIZE]: Old kind = 1 Old chunk size = 1
            [OpenMPTP][SCHEDULE_CHUNK_SIZE]: New kind = 1 New chunk size = 4""")
        self.assertIn(expected_output, self.output)

    def test_omp_numthreads_24_to_1(self):
        expected_output = textwrap.dedent("""\
            [OpenMPTP][NUMTHREADS]: Get old setting
            [OpenMPTP][NUMTHREADS]: Old Setting = 24
            [OpenMPTP][NUMTHREADS]: Set new setting
            [OpenMPTP][NUMTHREADS]: New Setting = 1""")
        self.assertIn(expected_output, self.output)

    def test_omp_shedule_type_1_to_1(self):
        expected_output = textwrap.dedent("""\
            [OpenMPTP][SCHEDULE_TYPE]: setting scheduling type
            [OpenMPTP][SCHEDULE_TYPE]: Old kind = 1 Old chunk size = 4
            [OpenMPTP][SCHEDULE_TYPE]: New kind = 1 New chunk size = 4""")
        self.assertIn(expected_output, self.output)

    def test_omp_chunk_size_4_to_4(self):
        expected_output = textwrap.dedent("""\
            [OpenMPTP][SCHEDULE_CHUNK_SIZE]: setting scheduling chunk size
            [OpenMPTP][SCHEDULE_CHUNK_SIZE]: Old kind = 1 Old chunk size = 4
            [OpenMPTP][SCHEDULE_CHUNK_SIZE]: New kind = 1 New chunk size = 4""")
        self.assertIn(expected_output, self.output)

    def test_omp_numthreads_1_to_1(self):
        expected_output = textwrap.dedent("""\
            [OpenMPTP][NUMTHREADS]: Get old setting
            [OpenMPTP][NUMTHREADS]: Old Setting = 1
            [OpenMPTP][NUMTHREADS]: Set new setting
            [OpenMPTP][NUMTHREADS]: New Setting = 1""")
        self.assertIn(expected_output, self.output)

    def test_no_exceptions(self):
        not_expected_output = textwrap.dedent("""\
            Uncaught exception with message:""")
        self.assertNotIn(not_expected_output, self.output)
        
        
        

if __name__ == '__main__':
    unittest.main()

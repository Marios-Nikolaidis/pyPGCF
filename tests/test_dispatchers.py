import unittest
from pypgcf import dispatchers

def _helper(x: int):
    return x+5

class TestModule(unittest.TestCase):
    def test_execute_command(self):
        cmd_pass = "pwd"
        cmd_fail = "LLL"
        ret_pass = dispatchers.execute_command(cmd_pass)
        ret_fail = dispatchers.execute_command(cmd_fail)

        self.assertEqual(ret_pass, 0)
        self.assertNotEqual(ret_fail, 0)

    def test_multiprocess_dispatch_callable_with_progress(self):
        inputs = list(range(5000))
        res = dispatchers.multiprocess_dispatch(_helper, inputs, 2, True, description="MP")
        # TODO: pytest doesn't show stdout
        total = sum(res)
        self.assertEqual(total, 12522500)

    def test_multiprocess_dispatch_callable_no_progress(self):
        inputs = list(range(5000))
        res = dispatchers.multiprocess_dispatch(_helper, inputs, 2, False, description="MP")
        # TODO: pytest doesn't show stdout
        total = sum(res)
        self.assertEqual(total, 12522500)

    def test_multiprocess_dispatch_system_with_progress(self):
        cmds = []
        for _ in range(100):
            cmds.append("echo TEST")
        res = dispatchers.multiprocess_dispatch("system", cmds, 2, True, description="Echo")
        total = len(res)
        self.assertEqual(total, 100)

    def test_multiprocess_dispatch_system_no_progress(self):
        cmds = []
        for _ in range(100):
            cmds.append("echo TEST")
        res = dispatchers.multiprocess_dispatch("system", cmds, 1, False, description="Echo")
        total = len(res)
        self.assertEqual(total, 100)

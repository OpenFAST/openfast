import os
from pathlib import Path
import shlex
import subprocess
import sys
import tempfile
import unittest

import openfastDrivers


class OpenfastDriversTest(unittest.TestCase):
    def setUp(self):
        self.temp_dir = tempfile.TemporaryDirectory()
        self.case_dir = Path(self.temp_dir.name) / "case"
        self.case_dir.mkdir()
        self.fixture = self.case_dir / "fixture.py"
        self.fixture.write_text(
            "import pathlib\n"
            "import os\n"
            "import sys\n"
            "print('stdout marker')\n"
            "print('stderr marker', file=sys.stderr)\n"
            "pathlib.Path(__file__).with_name('arguments.txt').write_text('\\n'.join(sys.argv[1:]))\n"
            "raise SystemExit(int(os.environ['OPENFAST_FIXTURE_EXIT_CODE']))\n"
        )

    def tearDown(self):
        self.temp_dir.cleanup()

    def invoke(self, *, exit_code=0, mode="normal", extra_flags=""):
        result = subprocess.run(
            [
                sys.executable,
                str(Path(__file__).resolve()),
                "--invoke",
                mode,
                str(exit_code),
                self.fixture.name,
                extra_flags,
            ],
            cwd=self.case_dir,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            check=False,
        )
        return result.returncode, result.stdout, result.stderr

    def test_normal_case_preserves_exit_code_flags_stdout_and_overwrite(self):
        for exit_code in (0, 7):
            with self.subTest(exit_code=exit_code):
                log_file = self.case_dir / "case.log"
                log_file.write_text("old log\n")

                return_code, output, errors = self.invoke(exit_code=exit_code, extra_flags="extra flag")

                self.assertEqual(return_code, exit_code)
                self.assertNotIn("stdout marker", output)
                self.assertNotIn("stdout marker", errors)
                log = log_file.read_text()
                self.assertIn("stdout marker", log)
                self.assertNotIn("old log", log)
                arguments = (self.case_dir / "arguments.txt").read_text().splitlines()
                self.assertEqual(arguments, ["input.fst", "extra", "flag"])

    def test_normal_case_captures_stderr_and_hides_child_streams(self):
        for exit_code in (0, 7):
            with self.subTest(exit_code=exit_code):
                log_file = self.case_dir / "case.log"
                log_file.write_text("old log\n")

                return_code, output, errors = self.invoke(exit_code=exit_code, extra_flags="extra flag")

                self.assertEqual(return_code, exit_code)
                self.assertIn("stderr marker", log_file.read_text())
                self.assertNotIn("stderr marker", output)
                self.assertNotIn("stderr marker", errors)

    def test_restart_case_preserves_exit_code_arguments_stdout_and_overwrite(self):
        for exit_code in (0, 7):
            with self.subTest(exit_code=exit_code):
                log_file = self.case_dir / "case_2.log"
                log_file.write_text("old log\n")
                return_code, output, errors = self.invoke(exit_code=exit_code, mode="restart")

                self.assertEqual(return_code, exit_code)
                self.assertNotIn("stdout marker", output)
                self.assertNotIn("stdout marker", errors)
                log = log_file.read_text()
                self.assertIn("stdout marker", log)
                self.assertNotIn("old log", log)
                self.assertEqual(
                    (self.case_dir / "arguments.txt").read_text().splitlines(),
                    ["-restart", "input"],
                )

    def test_restart_case_captures_stderr_and_hides_child_streams(self):
        for exit_code in (0, 7):
            with self.subTest(exit_code=exit_code):
                log_file = self.case_dir / "case_2.log"
                log_file.write_text("old log\n")
                return_code, output, errors = self.invoke(exit_code=exit_code, mode="restart")

                self.assertEqual(return_code, exit_code)
                self.assertIn("stderr marker", log_file.read_text())
                self.assertNotIn("stderr marker", output)
                self.assertNotIn("stderr marker", errors)

    def test_verbose_case_preserves_streams_arguments_and_log_state(self):
        log_file = self.case_dir / "case.log"
        for existing_log in (True, False):
            for exit_code in (0, 7):
                with self.subTest(existing_log=existing_log, exit_code=exit_code):
                    if existing_log:
                        log_file.write_text("old log\n")
                    elif log_file.exists():
                        log_file.unlink()

                    return_code, output, errors = self.invoke(exit_code=exit_code, mode="verbose", extra_flags="extra")

                    self.assertEqual(return_code, exit_code)
                    self.assertIn("stdout marker", output)
                    self.assertIn("stderr marker", errors)
                    self.assertNotIn("stderr marker", output)
                    self.assertNotIn("stdout marker", errors)
                    self.assertEqual(log_file.exists(), existing_log)
                    if existing_log:
                        self.assertEqual(log_file.read_text(), "old log\n")
                    self.assertEqual(
                        (self.case_dir / "arguments.txt").read_text().splitlines(),
                        ["input.fst", "extra"],
                    )


def invoke_fixture(mode, exit_code, fixture, extra_flags):
    os.environ["OPENFAST_FIXTURE_EXIT_CODE"] = str(exit_code)
    if os.name == "nt":
        executable = subprocess.list2cmdline([sys.executable, fixture])
    else:
        executable = " ".join((shlex.quote(sys.executable), shlex.quote(fixture)))
    log_file = {"normal": "case.log", "restart": "case_2.log"}.get(mode)
    return openfastDrivers._runCase(
        executable,
        "input.fst",
        log_file,
        sys.stdout,
        restart=mode in ("restart", "verbose"),
        ExtraFlags=extra_flags,
    )


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "--invoke":
        sys.exit(invoke_fixture(*sys.argv[2:]))
    unittest.main()

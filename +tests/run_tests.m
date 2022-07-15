%run unit tests and generate coverage report
import matlab.unittest.TestSuite
import matlab.unittest.TestRunner
import matlab.unittest.plugins.CodeCoveragePlugin

suite = TestSuite.fromPackage('tests');

runner = TestRunner.withTextOutput;

runner.addPlugin(CodeCoveragePlugin.forFolder('.'))
result = runner.run(suite);
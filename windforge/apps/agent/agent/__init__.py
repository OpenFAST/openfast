"""WindForge Agent -- local OpenFAST simulation executor.

The agent connects to the WindForge API server over WebSocket, receives
simulation jobs, executes them via the OpenFAST shared library, and streams
results back in real time.
"""

__version__ = "0.1.0"

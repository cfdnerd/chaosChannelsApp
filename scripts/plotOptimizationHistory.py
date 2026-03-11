#!/usr/bin/env python3
"""Plot optimization history from optimization.hst.

Supports:
- one-shot plotting after a run
- live plotting while optimization is running
"""

from __future__ import annotations

import argparse
import time
import textwrap
from pathlib import Path
from typing import List, Sequence, Tuple

EXCLUDED_PLOT_COLUMNS = {"PD0Sug"}


def parse_history(path: Path) -> Tuple[List[str], List[List[float]]]:
    """Parse optimization.hst, handling comments and repeated headers."""
    if not path.exists():
        return [], []

    header: List[str] = []
    rows_by_iter = {}

    with path.open("r", encoding="utf-8", errors="ignore") as f:
        for raw in f:
            line = raw.strip()
            if not line or line.startswith("#"):
                continue

            parts = line.split()
            if not parts:
                continue

            if parts[0] == "Iter":
                header = parts
                continue

            if not header:
                # Defer parsing until header appears.
                continue

            expected = len(header)
            if len(parts) < expected:
                # Likely an in-progress write; skip partial line.
                continue

            try:
                it = int(float(parts[0]))
                values = [float(v) for v in parts[1:expected]]
            except ValueError:
                continue

            rows_by_iter[it] = [float(it)] + values

    if not header or not rows_by_iter:
        return header, []

    ordered = [rows_by_iter[k] for k in sorted(rows_by_iter)]
    return header, ordered


def parse_column_selection(header: Sequence[str], columns_arg: str | None) -> List[int]:
    if not header:
        return []

    default_indices = [
        idx for idx in range(1, len(header)) if header[idx] not in EXCLUDED_PLOT_COLUMNS
    ]
    if not columns_arg:
        return default_indices

    wanted = [c.strip() for c in columns_arg.split(",") if c.strip()]
    if not wanted or wanted == ["all"]:
        return default_indices

    name_to_idx = {name: idx for idx, name in enumerate(header)}
    indices: List[int] = []
    for name in wanted:
        if name not in name_to_idx:
            raise ValueError(
                f"Unknown column '{name}'. Available: {', '.join(header[1:])}"
            )
        if name == "Iter":
            continue
        if name in EXCLUDED_PLOT_COLUMNS:
            continue
        indices.append(name_to_idx[name])

    if not indices:
        raise ValueError("No plottable columns selected.")

    return sorted(set(indices))


def build_figure(n_plots: int):
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(n_plots, 1, sharex=True, figsize=(11, 2.6 * n_plots))
    if n_plots == 1:
        axes = [axes]
    return fig, axes


def render(
    fig,
    axes,
    header: Sequence[str],
    rows: Sequence[Sequence[float]],
    col_indices: Sequence[int],
    src_path: Path,
) -> None:
    import matplotlib.pyplot as plt

    if not rows:
        return

    x = [r[0] for r in rows]

    for ax, idx in zip(axes, col_indices):
        y = [r[idx] for r in rows]
        ax.clear()
        ax.plot(x, y, "-", linewidth=1.6)
        ax.scatter([x[-1]], [y[-1]], s=18)
        ax.set_ylabel(header[idx])
        ax.grid(True, alpha=0.25)

    axes[-1].set_xlabel("Iteration")
    fig.suptitle(f"Optimization History: {src_path}")
    fig.tight_layout()
    plt.draw()


def main() -> int:
    epilog = textwrap.dedent(
        """\
        Examples:
          plotOptimizationHistory.py
          plotOptimizationHistory.py --live
          plotOptimizationHistory.py --live --interval 1.0
          plotOptimizationHistory.py -c PDCon,MeanT,MaxT
          plotOptimizationHistory.py --no-show --save solverLogs/optimization_plot.png
          plotOptimizationHistory.py -f solverLogs/optimization.hst --help-all
        """
    )
    parser = argparse.ArgumentParser(
        description="Plot optimization.hst parameters vs iteration.",
        epilog=epilog,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "-f",
        "--file",
        default="solverLogs/optimization.hst",
        help="Path to optimization.hst (default: solverLogs/optimization.hst)",
    )
    parser.add_argument(
        "-c",
        "--columns",
        default=None,
        help="Comma-separated columns to plot (default: all, excluding Iter and PD0Sug)",
    )
    parser.add_argument(
        "--live",
        action="store_true",
        help="Refresh plot continuously while the file is being updated.",
    )
    parser.add_argument(
        "--interval",
        type=float,
        default=2.0,
        help="Refresh interval in seconds for --live (default: 2.0)",
    )
    parser.add_argument(
        "--save",
        default=None,
        help="Optional output image path, e.g. solverLogs/optimization.png",
    )
    parser.add_argument(
        "--no-show",
        action="store_true",
        help="Run without opening a GUI window.",
    )
    parser.add_argument(
        "--help-all",
        action="store_true",
        help="Print extended help with workflow and troubleshooting notes.",
    )

    args = parser.parse_args()
    src = Path(args.file)

    if args.help_all:
        header, _ = parse_history(src)
        columns_msg = ", ".join(header[1:]) if header else "(none detected yet)"
        print(
            textwrap.dedent(
                f"""\
                Full Usage Guide
                ----------------
                Script:
                  scripts/plotOptimizationHistory.py

                Typical workflows:
                  1) Live during optimization:
                     ./scripts/plotOptimizationHistory.py --live
                  2) Post-run interactive:
                     ./scripts/plotOptimizationHistory.py
                  3) Post-run headless save:
                     ./scripts/plotOptimizationHistory.py --no-show --save solverLogs/optimization_plot.png

                Options:
                  -f, --file PATH        Input history file (default: solverLogs/optimization.hst)
                  -c, --columns LIST     Comma-separated columns (default: all except Iter and PD0Sug)
                  --live                 Refresh continuously until Ctrl+C
                  --interval SEC         Refresh interval for --live (default: 2.0)
                  --save PATH            Save plot image (png/jpg/pdf...)
                  --no-show              Do not open GUI window
                  -h, --help             Standard help
                  --help-all             This extended help

                Column selection:
                  - Use exact header names from optimization.hst.
                  - Example: --columns PDCon,MeanT,MaxT
                  - PD0Sug is intentionally excluded from plotting.
                  - Detected columns in {src}: {columns_msg}

                Notes:
                  - Repeated headers in optimization.hst are handled automatically.
                  - Incomplete trailing lines during active writes are skipped safely.
                  - With --save in --live mode, the output image is overwritten each refresh.
                """
            )
        )
        return 0

    if args.no_show:
        import matplotlib

        matplotlib.use("Agg")

    try:
        import matplotlib.pyplot as plt
    except Exception as exc:
        print(f"Failed to import matplotlib: {exc}")
        print("Install matplotlib, e.g. `pip install matplotlib`.")
        return 1

    header, rows = parse_history(src)
    if not header:
        print(f"No header found in {src}. Waiting for log content..." if args.live else f"No header found in {src}.")
        if not args.live:
            return 1

    try:
        col_indices = parse_column_selection(header, args.columns) if header else []
    except ValueError as exc:
        print(str(exc))
        return 1
    if header and not col_indices:
        print("No columns selected to plot.")
        return 1

    fig = None
    axes = None
    last_signature = None

    def update_once() -> bool:
        nonlocal fig, axes, last_signature, header, col_indices
        header_new, rows_new = parse_history(src)
        if not header_new or not rows_new:
            return False

        if header != header_new:
            header = header_new
            try:
                col_indices = parse_column_selection(header, args.columns)
            except ValueError as exc:
                print(str(exc))
                return False
            fig = None
            axes = None

        signature = (len(rows_new), int(rows_new[-1][0]))
        if signature == last_signature:
            return False
        last_signature = signature

        if fig is None:
            fig, axes = build_figure(len(col_indices))

        render(fig, axes, header, rows_new, col_indices, src)
        if args.save:
            Path(args.save).parent.mkdir(parents=True, exist_ok=True)
            fig.savefig(args.save, dpi=170, bbox_inches="tight")
        return True

    if args.live:
        print(f"Live plotting {src} every {args.interval:.2f}s (Ctrl+C to stop).")
        try:
            while True:
                update_once()
                if args.no_show:
                    time.sleep(max(args.interval, 0.2))
                else:
                    plt.pause(max(args.interval, 0.2))
        except KeyboardInterrupt:
            print("Stopped live plotting.")
            return 0

    changed = update_once()
    if not changed:
        print(f"No data rows parsed from {src}.")
        return 1

    if args.save:
        print(f"Saved plot: {args.save}")

    if not args.no_show:
        plt.show()

    return 0


if __name__ == "__main__":
    raise SystemExit(main())

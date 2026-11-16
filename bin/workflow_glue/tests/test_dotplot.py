"""Tests for PAF parsing, static PNG dotplot generation, and ezcharts dotplot."""
from pathlib import Path
import pytest


MOCK_PAF_DATA = (
    "hap1\t1000\t100\t900\t+\tchr16\t5000\t2000\t2800\t700\t800\t60\ttp:A:P\n"
    "hap2\t1000\t50\t950\t-\tchr16\t5000\t3000\t3900\t800\t900\t60\ttp:A:P\n"
)


def test_parse_paf_coordinates(tmp_path: Path):
    """Test parsing PAF alignments into coordinate segments for forward and reverse strands."""
    from workflow_glue.dotplot import parse_paf, compute_segments

    paf_file = tmp_path / "test.paf"
    paf_file.write_text(MOCK_PAF_DATA)

    records = parse_paf(paf_file)
    assert len(records) == 2

    # Record 1: hap1 on '+' strand
    r1 = records[0]
    assert r1["query_name"] == "hap1"
    assert r1["strand"] == "+"
    assert r1["target_name"] == "chr16"

    # Record 2: hap2 on '-' strand
    r2 = records[1]
    assert r2["query_name"] == "hap2"
    assert r2["strand"] == "-"

    segments = compute_segments(records)
    assert len(segments) == 2

    # Forward strand segment: x0=t_start, y0=q_start; x1=t_end, y1=q_end
    s1 = segments[0]
    assert s1["x0"] == 2000
    assert s1["y0"] == 100
    assert s1["x1"] == 2800
    assert s1["y1"] == 900
    assert s1["query_name"] == "hap1"

    # Reverse strand segment: x0=t_start, y0=q_end; x1=t_end, y1=q_start
    s2 = segments[1]
    assert s2["x0"] == 3000
    assert s2["y0"] == 950
    assert s2["x1"] == 3900
    assert s2["y1"] == 50
    assert s2["query_name"] == "hap2"


def test_generate_dotplot_png(tmp_path: Path):
    """Test static PNG dotplot generation using matplotlib."""
    from workflow_glue.dotplot import generate_dotplot_png

    paf_file = tmp_path / "test.paf"
    paf_file.write_text(MOCK_PAF_DATA)

    png_out = tmp_path / "test.dotplot.png"
    result_path = generate_dotplot_png(paf_file, png_out, alias="test_sample")

    assert result_path.exists()
    assert result_path.stat().st_size > 1000
    with open(result_path, "rb") as f:
        header = f.read(8)
    assert header == b"\x89PNG\r\n\x1a\n"


def test_create_ezcharts_dotplot(tmp_path: Path):
    """Test interactive ezcharts dotplot creation."""
    from workflow_glue.dotplot import create_ezcharts_dotplot
    from ezcharts.components.ezchart import EZChart

    paf_file = tmp_path / "test.paf"
    paf_file.write_text(MOCK_PAF_DATA)

    plot = create_ezcharts_dotplot(paf_file, alias="test_sample")
    assert plot is not None
    # Verify it can be wrapped into an EZChart component
    chart = EZChart(plot, "epi2melabs")
    html_str = str(chart)
    assert "bokeh" in html_str.lower()


def test_empty_paf_dotplot(tmp_path: Path):
    """Test handling of empty PAF files gracefully."""
    from workflow_glue.dotplot import parse_paf, generate_dotplot_png, create_ezcharts_dotplot

    empty_paf = tmp_path / "empty.paf"
    empty_paf.write_text("")

    records = parse_paf(empty_paf)
    assert records == []

    png_out = tmp_path / "empty.dotplot.png"
    result_path = generate_dotplot_png(empty_paf, png_out, alias="empty_sample")
    assert result_path.exists()
    assert result_path.stat().st_size > 500

    plot = create_ezcharts_dotplot(empty_paf, alias="empty_sample")
    assert plot is None

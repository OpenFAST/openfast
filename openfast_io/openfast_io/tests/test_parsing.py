"""Tests for openfast_io.parsing — especially fmt_field boundary cases."""
import pytest
from openfast_io.parsing import fmt_field, float_read, bool_read, int_read


class TestFmtField:
    """Boundary tests for fmt_field — the formatter that prevents overflow."""

    def test_short_string_pads_to_min_width(self):
        result = fmt_field(3.14)
        assert len(result) >= 28
        assert result.startswith('3.14')

    def test_long_value_extends_beyond_min_width(self):
        # A value whose str repr exceeds 28 chars
        long_val = -1.23456789012345678e-104
        result = fmt_field(long_val)
        assert len(result) >= len(str(long_val)) + 2

    def test_always_has_trailing_spaces(self):
        """Must always have at least 2 trailing spaces (prevents field concatenation)."""
        for val in [0, 1.0, -9999999999, 'default', 1e308, -1e-308, True, 'SomeFileName.dat']:
            result = fmt_field(val)
            assert result.endswith('  '), f'fmt_field({val!r}) lacks 2 trailing spaces: {result!r}'

    def test_integer_zero(self):
        result = fmt_field(0)
        assert result.strip() == '0'
        assert len(result) >= 28

    def test_large_positive_exponent(self):
        result = fmt_field(1.23e+200)
        assert '1.23e+200' in result or '1.23E+200' in result.upper()
        assert result.endswith('  ')

    def test_large_negative_exponent(self):
        result = fmt_field(-1.23e-200)
        assert result.rstrip() != result  # has trailing space
        assert '-1.23e-200' in result or '-1.23E-200' in result.upper()

    def test_string_default(self):
        result = fmt_field('default')
        assert result.startswith('default')
        assert len(result) >= 28

    def test_boolean_values(self):
        assert fmt_field(True).strip() == 'True'
        assert fmt_field(False).strip() == 'False'

    def test_none_value(self):
        result = fmt_field(None)
        assert result.strip() == 'None'

    def test_custom_min_width(self):
        result = fmt_field(42, min_width=10)
        assert len(result) >= 10

    def test_custom_min_width_smaller_than_value(self):
        result = fmt_field('a_very_long_filename_here.dat', min_width=5)
        # Should still have 2 trailing spaces
        assert result.endswith('  ')

    def test_negative_float_overflow_width_11(self):
        """The exact bug that affected AeroDyn blade tables: {:11} overflow."""
        val = -1.23456789e-04
        result = fmt_field(val, min_width=11)
        assert len(result) >= len(str(val)) + 2

    def test_stiffness_matrix_value(self):
        """The exact bug that affected BeamDyn: {:14} overflow for large stiffness."""
        val = 1.8634530e+10
        result = fmt_field(val, min_width=14)
        assert result.endswith('  ')
        assert len(result) >= len(str(val)) + 2


class TestFloatRead:
    def test_normal_float(self):
        assert float_read('3.14') == pytest.approx(3.14)

    def test_scientific_notation(self):
        assert float_read('1.5e-3') == pytest.approx(0.0015)

    def test_default_string(self):
        assert float_read('default') == 'default'
        assert float_read('DEFAULT') == 'DEFAULT'

    def test_trailing_comma(self):
        """The SubDyn GuyanDamp bug: trailing comma."""
        assert float_read('0.354293E+00') == pytest.approx(0.354293)

    def test_non_numeric(self):
        result = float_read('abc')
        assert result == 'abc'


class TestBoolRead:
    def test_true_variants(self):
        assert bool_read('True') is True
        assert bool_read('true') is True
        assert bool_read('T') is True
        assert bool_read('t') is True

    def test_false_variants(self):
        assert bool_read('False') is False
        assert bool_read('false') is False
        assert bool_read('F') is False

    def test_default(self):
        assert bool_read('default') == 'default'


class TestIntRead:
    def test_normal_int(self):
        assert int_read('42') == 42

    def test_default_string(self):
        assert int_read('default') == 'default'

    def test_non_numeric(self):
        assert int_read('abc') == 'abc'

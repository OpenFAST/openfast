import React, { useCallback, useState } from 'react';
import { Trash2, Plus, ChevronUp, ChevronDown } from 'lucide-react';
import clsx from 'clsx';

export interface DataGridColumn {
  name: string;
  key: string;
  type: 'number' | 'text' | 'select';
  min?: number;
  max?: number;
  step?: number;
  unit?: string;
  options?: { value: string; label: string }[];
  width?: string;
  precision?: number;
}

export interface DataGridProps<T extends Record<string, any>> {
  columns: DataGridColumn[];
  rows: T[];
  onChange: (rows: T[]) => void;
  onAddRow?: () => void;
  onDeleteRow?: (index: number) => void;
  rowHeight?: number;
  maxHeight?: string;
  readOnly?: boolean;
  className?: string;
}

export default function DataGrid<T extends Record<string, any>>({
  columns,
  rows,
  onChange,
  onAddRow,
  onDeleteRow,
  maxHeight = '500px',
  readOnly = false,
  className,
}: DataGridProps<T>) {
  const [selectedRows, setSelectedRows] = useState<Set<number>>(new Set());
  const [sortColumn, setSortColumn] = useState<string | null>(null);
  const [sortDirection, setSortDirection] = useState<'asc' | 'desc'>('asc');

  const handleCellChange = useCallback(
    (rowIndex: number, key: string, value: string | number) => {
      const updated = [...rows];
      const col = columns.find((c) => c.key === key);
      let parsed: any = value;

      if (col?.type === 'number') {
        parsed = value === '' ? 0 : parseFloat(value as string);
        if (isNaN(parsed)) parsed = 0;
        if (col.min !== undefined && parsed < col.min) parsed = col.min;
        if (col.max !== undefined && parsed > col.max) parsed = col.max;
      }

      updated[rowIndex] = { ...updated[rowIndex], [key]: parsed };
      onChange(updated);
    },
    [rows, columns, onChange]
  );

  const handleSelectRow = useCallback(
    (index: number, checked: boolean) => {
      const next = new Set(selectedRows);
      if (checked) {
        next.add(index);
      } else {
        next.delete(index);
      }
      setSelectedRows(next);
    },
    [selectedRows]
  );

  const handleSelectAll = useCallback(
    (checked: boolean) => {
      if (checked) {
        setSelectedRows(new Set(rows.map((_, i) => i)));
      } else {
        setSelectedRows(new Set());
      }
    },
    [rows]
  );

  const handleDeleteSelected = useCallback(() => {
    if (!onDeleteRow) return;
    const sortedIndices = Array.from(selectedRows).sort((a, b) => b - a);
    sortedIndices.forEach((i) => onDeleteRow(i));
    setSelectedRows(new Set());
  }, [selectedRows, onDeleteRow]);

  const handleSort = useCallback(
    (key: string) => {
      if (sortColumn === key) {
        setSortDirection((d) => (d === 'asc' ? 'desc' : 'asc'));
      } else {
        setSortColumn(key);
        setSortDirection('asc');
      }
    },
    [sortColumn]
  );

  const displayRows = React.useMemo(() => {
    if (!sortColumn) return rows.map((r, i) => ({ row: r, originalIndex: i }));
    const mapped = rows.map((r, i) => ({ row: r, originalIndex: i }));
    mapped.sort((a, b) => {
      const aVal = a.row[sortColumn] ?? 0;
      const bVal = b.row[sortColumn] ?? 0;
      const cmp = aVal < bVal ? -1 : aVal > bVal ? 1 : 0;
      return sortDirection === 'asc' ? cmp : -cmp;
    });
    return mapped;
  }, [rows, sortColumn, sortDirection]);

  const formatValue = (value: any, col: DataGridColumn): string => {
    if (value === null || value === undefined) return '';
    if (col.type === 'number' && col.precision !== undefined) {
      return Number(value).toFixed(col.precision);
    }
    return String(value);
  };

  return (
    <div className={clsx('flex flex-col', className)}>
      {/* Toolbar */}
      {!readOnly && (
        <div className="flex items-center gap-2 mb-2">
          {onAddRow && (
            <button
              onClick={onAddRow}
              className="flex items-center gap-1 px-3 py-1.5 text-xs font-medium
                         bg-blue-600 text-white rounded hover:bg-blue-700
                         transition-colors"
            >
              <Plus size={14} />
              Add Station
            </button>
          )}
          {onDeleteRow && selectedRows.size > 0 && (
            <button
              onClick={handleDeleteSelected}
              className="flex items-center gap-1 px-3 py-1.5 text-xs font-medium
                         bg-red-600 text-white rounded hover:bg-red-700
                         transition-colors"
            >
              <Trash2 size={14} />
              Delete ({selectedRows.size})
            </button>
          )}
          <span className="ml-auto text-xs text-gray-400">
            {rows.length} station{rows.length !== 1 ? 's' : ''}
          </span>
        </div>
      )}

      {/* Table */}
      <div
        className="overflow-auto border border-gray-700 rounded"
        style={{ maxHeight }}
      >
        <table className="w-full text-xs border-collapse">
          <thead className="sticky top-0 z-10">
            <tr className="bg-gray-800 text-gray-300">
              {!readOnly && onDeleteRow && (
                <th className="w-8 px-1 py-2 border-b border-gray-700">
                  <input
                    type="checkbox"
                    checked={
                      selectedRows.size === rows.length && rows.length > 0
                    }
                    onChange={(e) => handleSelectAll(e.target.checked)}
                    className="accent-blue-500"
                  />
                </th>
              )}
              <th className="w-10 px-2 py-2 border-b border-gray-700 text-left font-semibold text-gray-400">
                #
              </th>
              {columns.map((col) => (
                <th
                  key={col.key}
                  className="px-2 py-2 border-b border-gray-700 text-left font-semibold
                             cursor-pointer hover:bg-gray-750 select-none whitespace-nowrap"
                  style={{ width: col.width }}
                  onClick={() => handleSort(col.key)}
                >
                  <div className="flex items-center gap-1">
                    <span>{col.name}</span>
                    {col.unit && (
                      <span className="text-gray-500 font-normal">
                        [{col.unit}]
                      </span>
                    )}
                    {sortColumn === col.key && (
                      <span className="text-blue-400">
                        {sortDirection === 'asc' ? (
                          <ChevronUp size={12} />
                        ) : (
                          <ChevronDown size={12} />
                        )}
                      </span>
                    )}
                  </div>
                </th>
              ))}
            </tr>
          </thead>
          <tbody>
            {displayRows.length === 0 && (
              <tr>
                <td
                  colSpan={
                    columns.length +
                    1 +
                    (!readOnly && onDeleteRow ? 1 : 0)
                  }
                  className="px-4 py-8 text-center text-gray-500"
                >
                  No stations defined. Click "Add Station" to begin.
                </td>
              </tr>
            )}
            {displayRows.map(({ row, originalIndex }, displayIndex) => (
              <tr
                key={originalIndex}
                className={clsx(
                  'transition-colors',
                  displayIndex % 2 === 0 ? 'bg-gray-900' : 'bg-gray-850',
                  selectedRows.has(originalIndex) && 'bg-blue-900/30',
                  'hover:bg-gray-800'
                )}
              >
                {!readOnly && onDeleteRow && (
                  <td className="w-8 px-1 py-1 border-b border-gray-800">
                    <input
                      type="checkbox"
                      checked={selectedRows.has(originalIndex)}
                      onChange={(e) =>
                        handleSelectRow(originalIndex, e.target.checked)
                      }
                      className="accent-blue-500"
                    />
                  </td>
                )}
                <td className="w-10 px-2 py-1 border-b border-gray-800 text-gray-500 font-mono">
                  {originalIndex + 1}
                </td>
                {columns.map((col) => (
                  <td
                    key={col.key}
                    className="px-1 py-0.5 border-b border-gray-800"
                  >
                    {readOnly ? (
                      <span className="px-1 font-mono text-gray-300">
                        {formatValue(row[col.key], col)}
                      </span>
                    ) : col.type === 'select' ? (
                      <select
                        value={row[col.key] ?? ''}
                        onChange={(e) =>
                          handleCellChange(
                            originalIndex,
                            col.key,
                            e.target.value
                          )
                        }
                        className="w-full px-1 py-1 bg-gray-800 border border-gray-700
                                   rounded text-gray-200 text-xs focus:border-blue-500
                                   focus:outline-none focus:ring-1 focus:ring-blue-500/50"
                      >
                        <option value="">--</option>
                        {col.options?.map((opt) => (
                          <option key={opt.value} value={opt.value}>
                            {opt.label}
                          </option>
                        ))}
                      </select>
                    ) : (
                      <input
                        type={col.type === 'number' ? 'number' : 'text'}
                        value={
                          col.type === 'number'
                            ? row[col.key] ?? ''
                            : row[col.key] ?? ''
                        }
                        min={col.min}
                        max={col.max}
                        step={col.step ?? (col.type === 'number' ? 'any' : undefined)}
                        onChange={(e) =>
                          handleCellChange(
                            originalIndex,
                            col.key,
                            e.target.value
                          )
                        }
                        className="w-full px-1 py-1 bg-gray-800 border border-gray-700
                                   rounded text-gray-200 text-xs font-mono
                                   focus:border-blue-500 focus:outline-none
                                   focus:ring-1 focus:ring-blue-500/50
                                   [appearance:textfield]
                                   [&::-webkit-outer-spin-button]:appearance-none
                                   [&::-webkit-inner-spin-button]:appearance-none"
                      />
                    )}
                  </td>
                ))}
              </tr>
            ))}
          </tbody>
        </table>
      </div>
    </div>
  );
}

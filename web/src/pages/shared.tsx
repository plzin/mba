import * as React from 'react';
import { Dropdown } from '../ui/dropdown';
import Tooltip from '../ui/tooltip';
import { cn } from '../utils';
import { Formatter, SolutionAlgorithm } from '../mba';


interface InfoCardProps {
  children: React.ReactNode;
}

export function InfoCard({ children }: InfoCardProps) {
  return (
    <div className="card bg-info/30 dark:bg-info-content p-3 mb-4 text-sm text-info-content dark:text-info">{children}</div>
  );
}


interface OptionsLayoutProps {
  children: React.ReactNode;
  className?: string;
}

export function OptionsLayout({ children, className }: OptionsLayoutProps) {
  return (
    <div className={cn('grid px-2 grid-cols-1 sm:grid-cols-2 md:px-0 md:grid-cols-3 gap-4 pt-2', className)}>
      {children}
    </div>
  );
}

interface InputExpressionProps {
  expression: string;
  setExpression: (expression: string) => void;
  placeholder?: string;
  onEnter?: () => void;
}

export function InputExpression({
  expression,
  setExpression,
  onEnter,
  placeholder = 'Enter an expression',
}: InputExpressionProps) {
  return (
    <input
      type="text"
      placeholder={placeholder}
      className="input w-full"
      value={expression}
      onChange={e => setExpression(e.target.value)}
      onKeyDown={onEnter === undefined ? undefined : e => {
        if (e.key === 'Enter') {
          onEnter();
        }
      }}
    />
  );
}

interface IntegerWidthInputProps {
  numBits: number;
  setNumBits: (numBits: number) => void;
  tooltipSide?: 'top' | 'bottom' | 'left' | 'right';
}

export function IntegerWidthInput({
  numBits,
  setNumBits,
  tooltipSide = 'bottom',
}: IntegerWidthInputProps) {
  const [value, setValue] = React.useState(numBits.toString());
  return (
    <fieldset className="fieldset">
      <legend className="field-label">
        Integer width
        <Tooltip side={tooltipSide}>
          This can be any positive integer but for common integer widths (i.e. 8, 16, 32, 64, 128) the implementation is more efficient.
        </Tooltip>
      </legend>
      <input
        type="text"
        placeholder="Enter a number"
        className={cn('input', isNaN(Number(value)) && 'input-error')}
        value={value}
        onChange={(e) => {
          const value = e.target.value;
          setValue(value);
          if (!isNaN(Number(value))) {
            setNumBits(Number(value));
          }
        }}
      />
    </fieldset>
  );
}

interface OutputTypeDropdownProps {
  outputType: Formatter;
  setOutputType: (outputType: Formatter) => void;
}

export function OutputTypeDropdown({ outputType, setOutputType }: OutputTypeDropdownProps) {
  const outputTypes = [Formatter.C, Formatter.Rust, Formatter.Tex, Formatter.LLVM];
  return (
    <fieldset className="fieldset">
      <legend className="field-label">
        Output type
      </legend>
      <Dropdown
        entries={outputTypes}
        label={entry => Formatter[entry]}
        selected={outputType}
        setSelected={setOutputType}
        className="dropdown-center"
        buttonClassName="btn-secondary w-full"
        dropdownClassName="w-full"
      />
    </fieldset>
  );
}

interface SolutionAlgorithmDropdownProps {
  solutionAlgorithm: SolutionAlgorithm;
  setSolutionAlgorithm: (solutionAlgorithm: SolutionAlgorithm) => void;
}

export function SolutionAlgorithmDropdown({ solutionAlgorithm, setSolutionAlgorithm }: SolutionAlgorithmDropdownProps) {
  const solutionAlgorithms = [SolutionAlgorithm.Fast, SolutionAlgorithm.LeastComplexTerms, SolutionAlgorithm.ShortVector];
  const solutionAlgorithmLabels = {
    [SolutionAlgorithm.Fast]: 'Fast',
    [SolutionAlgorithm.LeastComplexTerms]: 'Least complex terms',
    [SolutionAlgorithm.ShortVector]: 'CVP',
  };
  return (
    <fieldset className="fieldset">
      <legend className="field-label">
        Solution algorithm
      </legend>
      <Dropdown
        entries={solutionAlgorithms}
        label={entry => solutionAlgorithmLabels[entry]}
        selected={solutionAlgorithm}
        setSelected={setSolutionAlgorithm}
        className="dropdown-center"
        buttonClassName="btn-secondary w-full"
        dropdownClassName="w-full"
      />
    </fieldset>
  );
}


interface ActionButtonProps {
  onClick: () => void;
  label: string;
  className?: string;
}
export function ActionButton({ onClick, label, className }: ActionButtonProps) {
  return (
    <fieldset className="fieldset">
      <button className={cn('btn btn-primary', className)} onClick={onClick}>{label}</button>
    </fieldset>
  );
}

interface CodeOutputProps {
  code?: string;
}
export function CodeOutput({ code }: CodeOutputProps) {
  if (!code) {
    return null;
  }
  return <div className="card mt-6 p-6 bg-base-200 overflow-x-scroll" dangerouslySetInnerHTML={{ __html: code }} />;
}

export interface PlaygroundButtonProps {
  label: string;
  onClick: () => void;
}

export function PlaygroundButton({ props }: { props?: PlaygroundButtonProps }) {
  if (!props) {
    return null;
  }
  return (
    <div className="flex flex-row mt-2">
      <button className="btn btn-warning ml-auto" onClick={props.onClick}>{props.label}</button>
    </div>
  );
}
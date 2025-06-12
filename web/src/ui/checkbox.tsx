import { cn } from '../utils';

interface LabeledCheckboxProps {
  children: React.ReactNode;
  checked: boolean;
  onChange: (checked: boolean) => void;
  className?: string;
  checkboxClassName?: string;
}

export function LabeledCheckbox({ children, checked, onChange, className, checkboxClassName }: LabeledCheckboxProps) {
  return (
    <label className={cn('label', className)}>
      <input type="checkbox" checked={checked} onChange={e => onChange(e.target.checked)} className={cn('checkbox', checkboxClassName)} />
      {children}
    </label>
  );
}
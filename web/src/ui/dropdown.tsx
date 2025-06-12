import { cn } from '../utils';

interface DropdownProps<T> {
  children?: React.ReactNode;
  entries: T[];
  label: (entry: T) => string;
  selected?: T;
  setSelected: (value: T) => void;
  className?: string;
  buttonClassName?: string;
  dropdownClassName?: string;
}

export function Dropdown<T>({ children, entries, label, selected, setSelected, className, buttonClassName, dropdownClassName }: DropdownProps<T>) {
  if (selected === undefined && children === undefined) {
    throw new Error('Either children or selected must be provided');
  }

  return (
    <details className={cn('dropdown', className)}>
      <summary className={cn('btn', buttonClassName)}>{children ?? label(selected!)}</summary>
      <ul className={cn('menu dropdown-content rounded-box bg-base-100 z-1 mt-3 p-2 shadow outline-1 outline-black/30', dropdownClassName)}>
        {entries.map((entry, index) => (
          <li key={index}>
            <a
              onClick={(e: React.MouseEvent<HTMLAnchorElement>) => {
                setSelected(entry);
                (e.target as HTMLAnchorElement).parentElement!.parentElement!.parentElement!.removeAttribute('open');
              }}
            >
              {label(entry)}
            </a>
          </li>
        ))}
      </ul>
    </details>
  );
}
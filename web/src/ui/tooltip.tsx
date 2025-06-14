import { cn } from '../utils';

interface TooltipProps {
  children: React.ReactNode;
  side?: 'top' | 'bottom' | 'left' | 'right';
}

export default function Tooltip({ children, side = 'bottom' }: TooltipProps) {
  const variants = {
    top: 'tooltip-top',
    bottom: 'tooltip-bottom',
    left: 'tooltip-left',
    right: 'tooltip-right',
  };

  return (
    <span className={cn('tooltip', variants[side])}>
      <div className="tooltip-content z-10">
        {children}
      </div>
      <svg xmlns="http://www.w3.org/2000/svg" width="14" height="14" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2" strokeLinecap="round" strokeLinejoin="round">
        <circle cx="12" cy="12" r="10"></circle>
        <path d="M12 16v-4"></path><path d="M12 8h.01"></path>
      </svg>
    </span>
  );
}
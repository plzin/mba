import React from 'react';
import { cn } from '../utils';

interface SliderProps {
  min: number;
  max: number;
  value: number;
  step?: number;
  onChange: (value: number) => void;
  className?: string;
  inputClassName?: string;
  showTicks?: 'all' | 'none' | 'bounds' | React.ReactNode[];
}

export function Slider({
  min,
  max,
  value,
  step = 1,
  onChange,
  className,
  inputClassName,
  showTicks = 'all',
}: SliderProps) {
  const handleChange = (event: React.ChangeEvent<HTMLInputElement>) => {
    onChange(Number(event.target.value));
  };

  const numSteps = Math.floor((max - min) / step) + 1;
  let ticks: React.ReactNode[] = [];
  if (showTicks === 'all') {
    ticks = Array.from({ length: numSteps }, (_, i) => (min + i * step).toString());
  } else if (showTicks === 'bounds') {
    ticks = [min.toString(), max.toString()];
  } else if (Array.isArray(showTicks)) {
    if (showTicks.length !== numSteps) {
      throw new Error(`Number of ticks must be ${numSteps}`);
    }
    ticks = showTicks;
  }

  return (
    <div className={className}>
      <input
        type="range"
        min={min}
        max={max}
        value={value}
        step={step}
        className={cn('range w-full', inputClassName)}
        onChange={handleChange}
      />
      {ticks.length > 0 && (
        <div className="flex justify-between px-2.5 mt-2 text-xs">
          {ticks.map((_, index) => (
            <span key={index}>|</span>
          ))}
        </div>
      )}
      {ticks.length > 0 && (
         <div className="flex justify-between px-2.5 mt-2 text-xs">
          {ticks.map((tick, index) => (
            <span key={index}>
              {tick}
            </span>
          ))}
        </div>
      )}
    </div>
  );
}

type LabeledSliderProps = SliderProps & {
  children: React.ReactNode;
};

export function LabeledSlider({
  children,
  ...props
}: LabeledSliderProps) {
  return (
    <fieldset className="fieldset">
      <legend className="field-label">
        {children}
      </legend>
      <Slider {...props} />
    </fieldset>
  );
}
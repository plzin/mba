import React from 'react';
import { InfoCard, InputExpression, IntegerWidthInput, OptionsLayout } from './shared';
import { fullLinearMBACheck, probabilisticLinearMBACheck } from '../mba';
import { cn } from '../utils';
import { LabeledSlider } from '../ui/slider';
import Tooltip from '../ui/tooltip';

export default function LinearChecker() {
  const [expression, setExpression] = React.useState('');
  const [numBits, setNumBits] = React.useState(8);
  const [numInputs, setNumInputs] = React.useState(0);
  const [error, setError] = React.useState<string>();
  const [isLinear, setIsLinear] = React.useState<boolean>();

  const updateExpression = (expression: string) => {
    setExpression(expression);
    setIsLinear(undefined);
    setError(undefined);
  };

  const check = () => {
    try {
      if (numInputs === 4) {
        setIsLinear(fullLinearMBACheck(expression, numBits));
      } else {
        setIsLinear(probabilisticLinearMBACheck(expression, Math.pow(10, numInputs + 2), numBits));
      }
    } catch (e) {
      if (typeof e === 'string') {
        setError(e);
      } else {
        setError('Unknown error');
        console.error(e);
      }
    }
  };

  const buttonColor = isLinear === undefined ? 'btn-primary' : isLinear ? 'btn-success' : 'btn-error';

  return (
    <div>
      <InfoCard>
        Determine (probabilistically) whether an expression can be represented as a linear MBA expression
        by checking if the fundamental theorem of mixed-boolean arithmetic holds for a number of random
        inputs. Note that setting "Number of inputs" to "All" is probably a bad idea.
      </InfoCard>
      <InputExpression expression={expression} setExpression={updateExpression} onEnter={check} />
      <OptionsLayout className="items-start">
        <IntegerWidthInput numBits={numBits} setNumBits={setNumBits} />
        <LabeledSlider
          min={0}
          max={4}
          value={numInputs}
          onChange={setNumInputs}
          showTicks={['100', '1000', '10000', '100000', <span className="text-error">All</span>]}
          inputClassName="range-secondary"
        >
          Number of inputs
          <Tooltip>
            The number of variable assignments to check. Setting this to "All" will check all possible assignments.
            You should probably not use this as even for 2 32-bit variables, this will have to check 2<sup>64</sup> assignments.
          </Tooltip>
        </LabeledSlider>
        <fieldset className="fieldset">
          <button className={cn('btn', buttonColor)} onClick={check}>Check</button>
        </fieldset>
      </OptionsLayout>
      {error && <p className="text-error">{error}</p>}
    </div>
  );
}
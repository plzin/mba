import * as React from 'react';
import { ActionButton, CodeOutput, InfoCard, InputExpression, IntegerWidthInput, OptionsLayout, OutputTypeDropdown, PlaygroundButton, SolutionAlgorithmDropdown, type PlaygroundButtonProps } from './shared';
import { deobfuscateLinear, Formatter, SolutionAlgorithm } from '../mba';
import { LabeledCheckbox } from '../ui/checkbox';
import Tooltip from '../ui/tooltip';
import { handleDisplayCode } from '../utils';

export default function DeobfuscateLinear() {
  const [expression, setExpression] = React.useState('');
  const [numBits, setNumBits] = React.useState(8);
  const [solutionAlgorithm, setSolutionAlgorithm] = React.useState<SolutionAlgorithm>(SolutionAlgorithm.LeastComplexTerms);
  const [detectBoolean, setDetectBoolean] = React.useState(false);
  const [outputType, setOutputType] = React.useState<Formatter>(Formatter.C);
  const [code, setCode] = React.useState<string>();
  const [playgroundButton, setPlaygroundButton] = React.useState<PlaygroundButtonProps>();

  const deobfuscate = () => {
    const getCode = () => deobfuscateLinear(
      expression, numBits, solutionAlgorithm, detectBoolean, outputType
    );
    handleDisplayCode(getCode, setCode, outputType, setPlaygroundButton);
  };

  return (
    <div>
      <InfoCard>
        Deobfuscates a linear MBA expression using the given algorithm.
      </InfoCard>
      <InputExpression
        expression={expression}
        setExpression={setExpression}
        onEnter={deobfuscate}
      />
      <OptionsLayout className="items-end">
        <IntegerWidthInput numBits={numBits} setNumBits={setNumBits} />
        <SolutionAlgorithmDropdown solutionAlgorithm={solutionAlgorithm} setSolutionAlgorithm={setSolutionAlgorithm} />
        <LabeledCheckbox
          checked={detectBoolean}
          onChange={setDetectBoolean}
          className="text-base-content text-sm font-medium self-center"
          checkboxClassName="checkbox-secondary"
        >
          Detect boolean
          <Tooltip>
            If enabled, the deobfuscation will detect if the expression can be expressed
            as a boolean expression and try to find a simple boolean expression for it without
            using the default deobfuscation algorithm that solves a linear system of equations.
          </Tooltip>
        </LabeledCheckbox>
        <OutputTypeDropdown outputType={outputType} setOutputType={setOutputType} />
        <ActionButton onClick={deobfuscate} label="Deobfuscate" />
      </OptionsLayout>
      <CodeOutput code={code} />
      <PlaygroundButton props={playgroundButton} />
    </div>
  );
}
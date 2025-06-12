import * as React from 'react';
import { ActionButton, CodeOutput, InfoCard, InputExpression, IntegerWidthInput, OptionsLayout, OutputTypeDropdown, PlaygroundButton, SolutionAlgorithmDropdown, type PlaygroundButtonProps } from './shared';
import { deobfuscate, Formatter, SolutionAlgorithm } from '../mba';
import { LabeledCheckbox } from '../ui/checkbox';
import Tooltip from '../ui/tooltip';
import { handleDisplayCode } from '../utils';

export default function Deobfuscate() {
  const [expression, setExpression] = React.useState('');
  const [numBits, setNumBits] = React.useState(8);
  const [solutionAlgorithm, setSolutionAlgorithm] = React.useState<SolutionAlgorithm>(SolutionAlgorithm.LeastComplexTerms);
  const [detectBoolean, setDetectBoolean] = React.useState(false);
  const [outputType, setOutputType] = React.useState<Formatter>(Formatter.C);
  const [code, setCode] = React.useState<string>();
  const [playgroundButton, setPlaygroundButton] = React.useState<PlaygroundButtonProps>();

  const onDeobfuscate = () => {
    const getCode = () => deobfuscate(
      expression, numBits, solutionAlgorithm, detectBoolean, outputType
    );
    handleDisplayCode(getCode, setCode, outputType, setPlaygroundButton);
  };

  return (
    <div>
      <InfoCard>
        Deobfuscates an expression by identifying the biggest linear MBA subexpressions and deobfuscating them individually.
        This can't handle non-linear MBA, e.g. the quadratic MBA obfuscation shown in the first blog post.
        As far as I know, the best way to deal with these is to try and invert the obfuscation, i.e. factor the expression.
        If you'd like to implement this, feel free to send a pull request!
      </InfoCard>
      <InputExpression
        expression={expression}
        setExpression={setExpression}
        onEnter={onDeobfuscate}
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
        <ActionButton onClick={onDeobfuscate} label="Deobfuscate" />
      </OptionsLayout>
      <CodeOutput code={code} />
      <PlaygroundButton props={playgroundButton} />
    </div>
  );
}
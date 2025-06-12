import React from 'react';
import { LabeledSlider } from '../ui/slider';
import { CodeOutput, InputExpression, IntegerWidthInput, ActionButton, OptionsLayout, OutputTypeDropdown, PlaygroundButton, type PlaygroundButtonProps, InfoCard } from './shared';
import { Formatter, obfuscate } from '../mba';
import { exampleShikiOptions, handleDisplayCode } from '../utils';

import { codeToHtml } from 'shiki';

const exampleExpr1 = await codeToHtml('x + y', exampleShikiOptions);
const exampleExpr2 = await codeToHtml('x & ~y', exampleShikiOptions);
const exampleExpr3 = await codeToHtml('1234', exampleShikiOptions);
const exampleExpr4 = await codeToHtml('x * y + z', exampleShikiOptions);

export default function Obfuscate() {
  const [expression, setExpression] = React.useState<string>('');
  const [outputType, setOutputType] = React.useState<Formatter>(Formatter.C);
  const [numBits, setNumBits] = React.useState(8);
  const [numAux, setNumAux] = React.useState(2);
  const [numRewrites, setNumRewrites] = React.useState(24);
  const [rewriteDepth, setRewriteDepth] = React.useState(3);
  const [code, setCode] = React.useState<string>();
  const [playgroundButton, setPlaygroundButton] = React.useState<PlaygroundButtonProps>();

  const onObfuscate = () => {
    const getCode = () => obfuscate(
      expression,
      numBits,
      numAux,
      rewriteDepth,
      numRewrites,
      outputType
    );
    handleDisplayCode(getCode, setCode, outputType, setPlaygroundButton);
  };

  return (
    <div>
      <InfoCard>
        Obfuscates an expression by identifying the biggest linear MBA subexpressions and rewriting them
        using linear MBA. This is just a simple example of how you might obfuscate an expression.
        Feel free to implement something more sophisticated and send a pull request!
      </InfoCard>
      <div className="tooltip tooltip-bottom w-full">
        <div className="tooltip-content z-50 max-w-128">
          This expression can be any mixed boolean-arithmetic expression, e.g.{' '}
          <span className="shiki" dangerouslySetInnerHTML={{ __html: exampleExpr1 }} />{', '}
          <span className="shiki" dangerouslySetInnerHTML={{ __html: exampleExpr2 }} />{', also constants '}
          <span className="shiki" dangerouslySetInnerHTML={{ __html: exampleExpr3 }} />.
          You can also use non-linear expressions such as{' '}
          <span className="shiki" dangerouslySetInnerHTML={{ __html: exampleExpr4 }} />, although the current
          implementation will just obfuscate the linear subexpressions and leave the non-linear parts as is.
        </div>
        <InputExpression expression={expression} setExpression={setExpression} onEnter={onObfuscate} />
      </div>
      <OptionsLayout className="items-end">
        <LabeledSlider
          min={0}
          max={12}
          value={numAux}
          onChange={setNumAux}
          inputClassName="range-secondary"
        >
          Auxiliary variables: {numAux}
        </LabeledSlider>
        <LabeledSlider
          min={4}
          max={1000}
          value={numRewrites}
          onChange={setNumRewrites}
          inputClassName="range-secondary"
          showTicks="bounds"
        >
          Rewrite operations: {numRewrites}
        </LabeledSlider>
        <LabeledSlider
          min={1}
          max={8}
          value={rewriteDepth}
          onChange={setRewriteDepth}
          inputClassName="range-secondary"
        >
          Rewrite operation depth: {rewriteDepth}
        </LabeledSlider>
        <IntegerWidthInput numBits={numBits} setNumBits={setNumBits} />
        <OutputTypeDropdown outputType={outputType} setOutputType={setOutputType} />
        <ActionButton onClick={onObfuscate} label="Obfuscate" />
      </OptionsLayout>
      <CodeOutput code={code} />
      <PlaygroundButton props={playgroundButton} />
    </div>
  );
}
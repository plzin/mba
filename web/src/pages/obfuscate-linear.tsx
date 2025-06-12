import React from 'react';
import { CodeOutput, InputExpression, IntegerWidthInput, ActionButton, OutputTypeDropdown, PlaygroundButton, type PlaygroundButtonProps, InfoCard } from './shared';
import Tooltip from '../ui/tooltip';
import { codeToHtml } from 'shiki';
import { Formatter, obfuscateLinear, parseAndPrintLBExpr } from '../mba';
import { handleDisplayCode } from '../utils';

interface RewriteOp {
  expr: string;
  html: string;
}

export default function ObfuscateLinear() {
  const [expression, setExpression] = React.useState<string>('');
  const [outputType, setOutputType] = React.useState<Formatter>(Formatter.C);
  const [numBits, setNumBits] = React.useState(8);
  const [rewriteOp, setRewriteOp] = React.useState<string>('');
  const [rewriteOpError, setRewriteOpError] = React.useState<string>();
  const [rewriteOps, setRewriteOps] = React.useState<RewriteOp[]>([]);
  const [code, setCode] = React.useState<string>();
  const [playgroundButton, setPlaygroundButton] = React.useState<PlaygroundButtonProps>();

  const addRewriteOperation = () => {
    let expr: string;
    try {
      expr = parseAndPrintLBExpr(rewriteOp);
    } catch (e) {
      if (typeof e === 'string') {
        setRewriteOpError(e);
      } else {
        setRewriteOpError('Unknown error while parsing');
      }
      return;
    }

    // We always highlight as C because that's what the expressions are printed as.
    // Even if they weren't, C supports both `~` and `!` unary operators, which is
    // the only difference for these kinds of expressions.
    codeToHtml(expr, {
      lang: 'c',
      themes: {
        light: 'github-light-default',
        dark: 'github-dark-default',
      },
      defaultColor: false,
      structure: 'inline',
    }).then(html => {
      setRewriteOps([...rewriteOps, { expr, html }]);
      setRewriteOpError(undefined);
    }).catch(e => {
      setRewriteOpError('Error while highlighting');
      console.error(e);
    });
  };

  const removeRewriteOperation = (index: number) => {
    setRewriteOps(rewriteOps.filter((_, i) => i !== index));
  };

  const obfuscate = () => {
    const getCode = () => obfuscateLinear(
      expression,
      rewriteOps.map(op => op.expr),
      numBits,
      true,
      outputType
    );
    handleDisplayCode(getCode, setCode, outputType, setPlaygroundButton);
  };

  return (
    <div>
      <InfoCard>
        Rewrites a linear mixed-boolean expression using the set of rewrite operations.
      </InfoCard>
      <InputExpression expression={expression} setExpression={setExpression} onEnter={obfuscate} />
      <div className="grid px-2 grid-cols-1 sm:grid-cols-[minmax(0,2fr)_minmax(0,1fr)] gap-4 pt-2">
        <div>
          <fieldset className="fieldset bg-base-200 border-base-300 rounded-box border p-4 pt-0">
            <legend className="fieldset-legend">
              Rewrite operations
              <Tooltip>
                <p>
                  The expressions that appear in the linear combination. They have to be linear MBA expressions themselves.
                </p>
              </Tooltip>
            </legend>
            <ul className="list">
              {rewriteOps.map((operation, index) => (
                <li className="list-row" key={index}>
                  <div className="list-col-grow flex flex-row gap-2 items-center">
                    <span className="shiki flex-1" dangerouslySetInnerHTML={{ __html: operation.html }} />
                    <button
                      className="btn btn-square btn-dash btn-error btn-xs"
                      onClick={() => removeRewriteOperation(index)}
                    >
                      x
                    </button>
                  </div>
                </li>
              ))}
              <li className="list-row">
                <div className="list-col-grow">
                  <div className="flex flex-row gap-2">
                    <input
                      type="text"
                      className="input flex-1"
                      placeholder="Rewrite operation"
                      value={rewriteOp}
                      onChange={e => setRewriteOp(e.target.value)}
                      onKeyDown={e => {
                        if (e.key === 'Enter') {
                          addRewriteOperation();
                        }
                      }}
                    />
                    <button className="btn btn-accent" onClick={addRewriteOperation}>Add</button>
                  </div>
                  {rewriteOpError && <p className="text-error">{rewriteOpError}</p>}
                </div>
              </li>
            </ul>
          </fieldset>
        </div>
        <div className="space-y-4 pt-4">
          <IntegerWidthInput numBits={numBits} setNumBits={setNumBits} />
          <OutputTypeDropdown outputType={outputType} setOutputType={setOutputType} />
          <ActionButton onClick={obfuscate} label="Obfuscate" />
        </div>
      </div>
      <CodeOutput code={code} />
      <PlaygroundButton props={playgroundButton} />
    </div>
  );
}
import * as React from 'react';
import { ActionButton, InfoCard, IntegerWidthInput } from './shared';
import { solveLinearSystem } from '../mba';
import { tex2html } from '../utils';

interface TexTrace {
  diag: string;
  scalarSystem: string;
  linearSolutions: string;
  vectorSolution: string;
  finalSolution: string;
}

export default function LinearSystem() {
  const [system, setSystem] = React.useState<string>('');
  const [numBits, setNumBits] = React.useState(8);
  const [trace, setTrace] = React.useState<TexTrace>();
  const [error, setError] = React.useState<string>();

  const systemHtml = React.useMemo(() => {
    let tex = '\\begin{aligned}';

    for (let line of system.split('\n')) {
      line = line.trim();
      if (line.length === 0) {
        continue;
      }

      const elements = line.split(' ');
      const rhs = elements.length > 1 ? elements.pop() : '';
      for (const [i, element] of elements.entries()) {
        tex += `${element}x_{${i+1}}&+`;
      }
      tex = tex.slice(0, -1);
      tex += `\\equiv ${rhs}&\\pmod{2^{${numBits}}}\\\\`;
    }
    tex += '\\end{aligned}';

    return tex2html(tex);
  }, [system, numBits]);

  const onSolve = () => {
    try {
       const trace = solveLinearSystem(system, numBits);
       setTrace({
        diag: tex2html(trace.diag),
        scalarSystem: tex2html(trace.scalarSystem),
        linearSolutions: tex2html(trace.linearSolutions),
        vectorSolution: tex2html(trace.vectorSolution),
        finalSolution: tex2html(trace.finalSolution),
       });
    } catch (e) {
      if (typeof e === 'string') {
        setError(e);
      } else {
        console.error(e);
        setError('Unknown error while solving the linear system. See console for details.');
      }
    }
  };

  return (
    <div>
      <InfoCard>
        <span>Solves a system of linear equations mod 2<sup>n</sup>.</span>
      </InfoCard>
      <div className="flex flex-row gap-4">
        <textarea
          className="textarea w-full"
          value={system}
          onChange={(e) => setSystem(e.target.value)}
          placeholder={'Enter the components of the linear system, e.g.\n3 4 5\n2 5 7'}
        />
        <div>
          <IntegerWidthInput numBits={numBits} setNumBits={setNumBits} tooltipSide="left" />
          <ActionButton onClick={onSolve} label="Solve" />
        </div>
      </div>
      <TexOutput tex={systemHtml} />
      {error && <div className="text-error">{error}</div>}
      {trace && <div>
        The diagonalization of the matrix is
        <TexOutput tex={trace.diag} />
        The new system is
        <TexOutput tex={trace.scalarSystem} />
        which results in the single variable linear congruences
        <TexOutput tex={trace.linearSolutions} />
        {trace.vectorSolution.length === 0 ? 'So overall the system has no solution.' :
        <>
        The vector form is
        <TexOutput tex={trace.vectorSolution} />
        The final solution to the original system is
        <TexOutput tex={trace.finalSolution} />
        </>
        }
      </div>}
    </div>
  );
}

function TexOutput({ tex }: { tex: string }) {
  return <div className="overflow-x-auto" dangerouslySetInnerHTML={{ __html: tex }} />;
}
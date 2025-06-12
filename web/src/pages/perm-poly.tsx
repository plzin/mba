import React from 'react';
import { ActionButton, InfoCard, InputExpression, IntegerWidthInput, OptionsLayout } from './shared';
import { invertPermutationPolynomial, randomPermutationPolynomial } from '../mba';
import { tex2html } from '../utils';

export default function PermPoly() {
  const [expression, setExpression] = React.useState('');
  const [numBits, setNumBits] = React.useState(8);
  const [html, setHtml] = React.useState('');

  const randomPolynomial = () => {
    setExpression(randomPermutationPolynomial(numBits));
  };

  const computeInverse = () => {
    try {
      setHtml(tex2html(invertPermutationPolynomial(expression, numBits)));
    } catch (e) {
      if (typeof e === 'string') {
        setHtml(`<span class="text-error">${e}</span>`);
      } else {
        setHtml('<span class="text-error">Unknown error. Check console.</span>');
        console.error(e);
      }
    }
  };

  return (
    <div>
      <InfoCard>
        Generate and invert binary permutation polynomials.
      </InfoCard>
      <InputExpression
        expression={expression}
        setExpression={setExpression}
        placeholder="Enter a permutation polynomial"
        onEnter={computeInverse}
      />
      <OptionsLayout className="items-end">
        <IntegerWidthInput numBits={numBits} setNumBits={setNumBits} />
        <ActionButton className="btn-secondary" label="Random" onClick={randomPolynomial} />
        <ActionButton label="Compute inverse" onClick={computeInverse} />
      </OptionsLayout>
      <div className="katex-wrap" dangerouslySetInnerHTML={{ __html: html }} />
    </div>
  );
}
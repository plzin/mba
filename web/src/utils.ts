import { codeToHtml } from 'shiki';
import type { PlaygroundButtonProps } from './pages/shared';
import { Formatter } from './mba';
import katex from 'katex';

/**
 * Concatenates class names into a single string.
 * This currently doesn't use twMerge, because I didn't want to add a dependency.
 * So it has to be used with care, but since this is a small project, it should be okay.
 * @param classes - The class names to concatenate.
 * @returns The concatenated class names.
 */
export function cn(...classes: (string | undefined | false)[]) {
  return classes.filter(Boolean).join(' ');
}

/**
 * Converts a LaTeX string to HTML.
 * @param tex - The LaTeX string to convert.
 * @returns The HTML string.
 */
export function tex2html(tex: string): string {
  if (tex.length === 0) {
    return '';
  }

  const html = katex.renderToString(tex, {
    displayMode: true,
    output: 'html',
    throwOnError: false,
  });

  return html;
}


/**
 * Handles the display of the code.
 * @param genCode - The function that generates the code.
 * @param setCode - The function that sets the code.
 * @param outputType - The output type.
 * @param setPlaygroundButton - The function that sets the playground button.
 */
export function handleDisplayCode(
  genCode: () => string,
  setCode: (html: string) => void,
  outputType: Formatter,
  setPlaygroundButton: (button: PlaygroundButtonProps | undefined) => void,
) {
  // Try to generate the code.
  let code: string;
  try {
    code = genCode();
  } catch (e) {
    if (typeof e === 'string') {
      setCode(`<span class="text-error">${e}</span>`);
    } else {
      console.error(e);
      setCode('<span class="text-error">Unknown error while obfuscationg. See console for details.</span>');
      setPlaygroundButton(undefined);
    }
    return;
  }

  // Post-process the code.
  code = postProcessCode(code);

  let playgroundButton: PlaygroundButtonProps | undefined = undefined;

  // Set the playground button.
  if (outputType == Formatter.C) {
    const args = code.split(',').map(() => '0').join(', ');
    const ceCode = encodeURIComponent(`#include <cstdint>\n#include <iostream>\n\n${code}\n\nint main() {\n\tstd::cout << +f(${args}) << "\\n";\n}`);
    playgroundButton = {
      label: 'Open in Compiler Explorer',
      onClick: () => {
        window.open(`https://godbolt.org/#g:!((g:!((g:!((h:codeEditor,i:(filename:'1',fontScale:14,fontUsePx:'0',j:1,lang:c%2B%2B,selection:(endColumn:1,endLineNumber:1,positionColumn:1,positionLineNumber:1,selectionStartColumn:1,selectionStartLineNumber:1,startColumn:1,startLineNumber:1),source:'${ceCode}'),l:'5',n:'0',o:'C%2B%2B+source+%231',t:'0')),k:50,l:'4',n:'0',o:'',s:0,t:'0'),(g:!((h:executor,i:(argsPanelShown:'1',compilationPanelShown:'0',compiler:clang_trunk,compilerOutShown:'0',execArgs:'',execStdin:'',fontScale:14,fontUsePx:'0',j:1,lang:c%2B%2B,libs:!(),options:'',source:1,stdinPanelShown:'1',tree:'1',wrap:'1'),l:'5',n:'0',o:'Executor+x86-64+clang+(trunk)+(C%2B%2B,+Editor+%231)',t:'0')),k:50,l:'4',n:'0',o:'',s:0,t:'0')),l:'2',n:'0',o:'',t:'0')),version:4`);
      },
    };
  } else if (outputType == Formatter.Rust) {
    const args = code.split(',').map(() => 'Wrapping(0)').join(', ');
    const pgCode = encodeURIComponent(`use std::num::Wrapping;\n\nfn main() {\n\tprintln!("{}", f(${args}));\n}\n\n${code}`);
    playgroundButton = {
      label: 'Open in Rust Playground',
      onClick: () => {
        window.open(`https://play.rust-lang.org/?version=stable&mode=release&edition=2024&code=${pgCode}`);
      },
    };
  }

  // Highlight the code.
  codeToHtml(code, {
    lang: Formatter[outputType].toLowerCase(),
    themes: {
      light: 'github-light-default',
      dark: 'github-dark-default',
    },
    defaultColor: false,
  }).then(html => {
    setCode(html);
    setPlaygroundButton(playgroundButton);
  }).catch(e => {
    console.error(e);
    setCode('<span class="text-error">Error while highlighting</span>');
    setPlaygroundButton(undefined);
  });
}

/**
 * Split the code into multiple lines, so that it's more readable.
 * @param code - The code to post-process.
 * @returns The post-processed code.
 */
function postProcessCode(code: string) {
  let s = '';
  const lines = code.split('\n');
  for (let l of lines) {
    let tabs = 0;
    while (tabs < l.length && l[tabs] == '\t') {
      tabs++;
    }
    l = l.substring(tabs);

    inner: for (let j = 0; ; j++) {
      let max = 68 - 4 * tabs;
      if (j == 0) {
        s += '\t'.repeat(tabs);
      } else {
        s += '\t'.repeat(tabs + 1);
        max -= 4;
      }

      if (l.length >= max) {
        for (let i = max; i >= 0; i--) {
          if (l[i] == ' ') {
            s += l.substring(0, i);
            s += '\n';
            l = l.substring(i + 1);
            continue inner;
          }
        }
      }

      // If we never found a space then just put the whole string into the line anyways.
      s += l;
      s += '\n';
      break;
    }
  }

  return s;
}
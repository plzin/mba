import React from 'react';
import Obfuscate from './pages/obfuscate';
import ObfuscateLinear from './pages/obfuscate-linear';
import Deobfuscate from './pages/deobfuscate';
import DeobfuscateLinear from './pages/deobfuscate-linear';
import LinearSystem from './pages/linear-system';
import ThemeController from './ui/theme-controller';
import { Dropdown } from './ui/dropdown';
import PermPoly from './pages/perm-poly';
import LinearChecker from './pages/linear-checker';

function App() {
  const pages = {
    'obfuscate': {
      title: 'Obfuscate',
      component: <Obfuscate />,
    },
    'obfuscate-linear': {
      title: 'Obfuscate Linear',
      component: <ObfuscateLinear />,
    },
    'deobfuscate': {
      title: 'Deobfuscate',
      component: <Deobfuscate />,
    },
    'deobfuscate-linear': {
      title: 'Deobfuscate Linear',
      component: <DeobfuscateLinear />,
    },
    'linear-checker': {
      title: 'Linear Checker',
      component: <LinearChecker />,
    },
    'linear-system': {
      title: 'Linear System',
      component: <LinearSystem />,
    },
    'perm-poly': {
      title: 'Binary Permutation Polynomials',
      component: <PermPoly />,
    }
  };

  const [currentPage, setCurrentPage] = React.useState(() => {
    // Get initial page from URL parameters
    const params = new URLSearchParams(window.location.search);
    return params.get('page') ?? 'obfuscate';
  });

  const navigateToPage = (page: string) => {
    const params = new URLSearchParams(window.location.search);
    params.set('page', page);
    window.history.pushState({}, '', `?${params.toString()}`);
    setCurrentPage(page);
  };

  // Listen to browser back/forward buttons
  React.useEffect(() => {
    const handlePopState = () => {
      const params = new URLSearchParams(window.location.search);
      setCurrentPage(params.get('page') ?? 'home');
    };

    window.addEventListener('popstate', handlePopState);
    return () => window.removeEventListener('popstate', handlePopState);
  }, []);

  const page = pages[currentPage as keyof typeof pages];

  return (
    <div className="sm:max-w-screen-md mx-auto font-mono">
      <div className="navbar">
        <div className="navbar-start w-full sm:w-1/2">
          <Dropdown
            entries={Object.entries(pages).map(([key, value]) => ({ key, value }))}
            label={entry => entry.value.title}
            setSelected={page => navigateToPage(page.key)}
            buttonClassName="btn-ghost btn-circle"
            dropdownClassName="w-56"
          >
            <svg xmlns="http://www.w3.org/2000/svg" className="h-5 w-5" fill="none" viewBox="0 0 24 24" stroke="currentColor">
              <path strokeLinecap="round" strokeLinejoin="round" strokeWidth="2" d="M4 6h16M4 12h16M4 18h7" />
            </svg>
          </Dropdown>
          <div className="text-xs sm:hidden">
            Mixed Boolean-Arithmetic
          </div>
        </div>
        <div className="navbar-center hidden sm:block sm:text-xl">
          Mixed Boolean-Arithmetic
        </div>
        <div className="navbar-end gap-3">
          <a href="docs/mba/index.html" title="Crate documentation">
            <svg xmlns="http://www.w3.org/2000/svg" className="h-5 w-5" fill="currentColor" viewBox="0 0 576 512">
              <path d="M290.8 48.6l78.4 29.7L288 109.5 206.8 78.3l78.4-29.7c1.8-.7 3.8-.7 5.7 0zM136 92.5l0 112.2c-1.3 .4-2.6 .8-3.9 1.3l-96 36.4C14.4 250.6 0 271.5 0 294.7L0 413.9c0 22.2 13.1 42.3 33.5 51.3l96 42.2c14.4 6.3 30.7 6.3 45.1 0L288 457.5l113.5 49.9c14.4 6.3 30.7 6.3 45.1 0l96-42.2c20.3-8.9 33.5-29.1 33.5-51.3l0-119.1c0-23.3-14.4-44.1-36.1-52.4l-96-36.4c-1.3-.5-2.6-.9-3.9-1.3l0-112.2c0-23.3-14.4-44.1-36.1-52.4l-96-36.4c-12.8-4.8-26.9-4.8-39.7 0l-96 36.4C150.4 48.4 136 69.3 136 92.5zM392 210.6l-82.4 31.2 0-89.2L392 121l0 89.6zM154.8 250.9l78.4 29.7L152 311.7 70.8 280.6l78.4-29.7c1.8-.7 3.8-.7 5.7 0zm18.8 204.4l0-100.5L256 323.2l0 95.9-82.4 36.2zM421.2 250.9c1.8-.7 3.8-.7 5.7 0l78.4 29.7L424 311.7l-81.2-31.1 78.4-29.7zM523.2 421.2l-77.6 34.1 0-100.5L528 323.2l0 90.7c0 3.2-1.9 6-4.8 7.3z"/>
            </svg>
          </a>
          <a href="https://github.com/plzin/mba" className="-mt-0.5" title="GitHub">
            <svg xmlns="http://www.w3.org/2000/svg" className="h-5 w-5" fill="currentColor" viewBox="0 0 16 16">
              <path d="M8 0C3.58 0 0 3.58 0 8c0 3.54 2.29 6.53 5.47 7.59.4.07.55-.17.55-.38 0-.19-.01-.82-.01-1.49-2.01.37-2.53-.49-2.69-.94-.09-.23-.48-.94-.82-1.13-.28-.15-.68-.52-.01-.53.63-.01 1.08.58 1.23.82.72 1.21 1.87.87 2.33.66.07-.52.28-.87.51-1.07-1.78-.2-3.64-.89-3.64-3.95 0-.87.31-1.59.82-2.15-.08-.2-.36-1.02.08-2.12 0 0 .67-.21 2.2.82.64-.18 1.32-.27 2-.27s1.36.09 2 .27c1.53-1.04 2.2-.82 2.2-.82.44 1.1.16 1.92.08 2.12.51.56.82 1.27.82 2.15 0 3.07-1.87 3.75-3.65 3.95.29.25.54.73.54 1.48 0 1.07-.01 1.93-.01 2.2 0 .21.15.46.55.38A8.01 8.01 0 0 0 16 8c0-4.42-3.58-8-8-8" />
            </svg>
          </a>
          <ThemeController />
        </div>
      </div>
      <div className="p-2">
        {page ? page.component : <div>Page not found</div>}
      </div>
    </div>
  );
}

export default App;

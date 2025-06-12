// Initialize the wasm module and export the functions that are defined in the wasm module.

import init from './wasm';

const wasm = await init();
wasm.setPanicHook();

export * from './wasm';
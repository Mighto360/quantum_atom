build:
	cargo build --target wasm32-unknown-unknown

start: build
	xcopy .\target\wasm32-unknown-unknown\debug\quantum_atom.wasm .\www\bin\ /y
	http-server .\www -c-1
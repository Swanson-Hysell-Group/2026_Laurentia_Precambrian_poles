.PHONY: help install build build-execute start clean api-docs

BUILD_DIR = _build

help: ## Show this help
	@grep -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) | \
		awk 'BEGIN {FS = ":.*?## "}; {printf "  \033[36m%-20s\033[0m %s\n", $$1, $$2}'

install: ## Create mamba environment from environment.yml
	mamba env create -f environment.yml --yes
	@echo ""
	@echo "Environment created. Activate with:"
	@echo "  mamba activate laurentia-poles"

build: ## Build the book (static HTML)
	jupyter book build --html

build-execute: ## Build the book with notebook execution
	jupyter book build --html --execute

start: ## Start local dev server with live reload
	myst start

clean: ## Remove build artifacts
	jupyter book clean
	rm -rf $(BUILD_DIR)

api-docs: ## Regenerate pole_tools API reference
	python scripts/generate_api_docs.py

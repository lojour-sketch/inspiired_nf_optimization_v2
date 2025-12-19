#!/usr/bin/env python3
import yaml
import sys
import json

if len(sys.argv) != 2:
    print("Usage: python debug_yaml.py <meta.yml>")
    sys.exit(1)

with open(sys.argv[1], 'r') as f:
    data = yaml.safe_load(f)

print("=== YAML Structure ===")
print(f"Keys found: {list(data.keys())}")
print()

if 'output' in data:
    print("=== OUTPUT SECTION ===")
    print(json.dumps(data['output'], indent=2, default=str))
else:
    print("❌ NO OUTPUT SECTION FOUND")
    
print("\n=== Full YAML ===")
print(json.dumps(data, indent=2, default=str))
#!/usr/bin/env python3
"""Quick script to get precision scores for all classes."""

import subprocess
import re

# Get list of classes
result = subprocess.run(['uv', 'run', 'autarch', 'list-classes'], 
                       capture_output=True, text=True)

classes = []
for line in result.stdout.split('\n'):
    if '│' in line and 'autarch.ontology' in line:
        class_name = line.split('│')[1].strip()
        if class_name and class_name != 'Class Name':
            classes.append(class_name)

print("Class\tPrecision")
precision_scores = []

for class_name in classes:
    try:
        result = subprocess.run(['uv', 'run', 'autarch', 'eval', class_name], 
                              capture_output=True, text=True, timeout=30)
        
        # Extract precision from output
        precision_match = re.search(r'│ Precision\s+│\s+([\d.]+)', result.stdout)
        if precision_match:
            precision = float(precision_match.group(1))
            precision_scores.append((class_name, precision))
            print(f"{class_name}\t{precision:.3f}")
        else:
            print(f"{class_name}\tError")
    except subprocess.TimeoutExpired:
        print(f"{class_name}\tTimeout")

print("\nWorst performers (lowest precision):")
precision_scores.sort(key=lambda x: x[1])
for class_name, precision in precision_scores[:10]:
    print(f"{class_name}: {precision:.3f}")

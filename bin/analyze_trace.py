#!/usr/bin/env python3
import sys
import csv
import re
from datetime import timedelta

def parse_time(time_str):
    """Convert Nextflow time format to seconds
    Handles formats like: '1m 30s', '1m30s', '1m', '30s', '1.5h', '100ms', '-'
    """
    if not time_str or time_str == '-' or time_str == '':
        return 0
    
    time_str = time_str.strip()
    total_seconds = 0
    
    # Handle combined formats like "1m 30s" or "1h 2m 3s"
    # Find all number+unit pairs
    pattern = r'([\d.]+)\s*([a-z]+)'
    matches = re.findall(pattern, time_str.lower())
    
    if not matches:
        return 0
    
    conversions = {
        'ms': 0.001,
        's': 1,
        'm': 60,
        'h': 3600,
        'd': 86400
    }
    
    for value_str, unit in matches:
        try:
            value = float(value_str)
            multiplier = conversions.get(unit, 0)
            total_seconds += value * multiplier
        except ValueError:
            continue
    
    return total_seconds

def format_time(seconds):
    """Format seconds into a readable string"""
    if seconds < 60:
        return f"{seconds:.2f}s"
    elif seconds < 3600:
        return f"{seconds/60:.2f}m ({seconds:.0f}s)"
    elif seconds < 86400:
        hours = seconds / 3600
        return f"{hours:.2f}h ({seconds/60:.1f}m)"
    else:
        days = seconds / 86400
        return f"{days:.2f}d ({seconds/3600:.1f}h)"

def analyze_trace(trace_file):
    executed_time = 0
    cached_time = 0
    executed_count = 0
    cached_count = 0
    
    try:
        with open(trace_file, 'r') as f:
            reader = csv.DictReader(f, delimiter='\t')
            
            for row in reader:
                realtime = parse_time(row.get('realtime', '0'))
                is_cached = row.get('cached', '0') == '1'
                
                if is_cached:
                    cached_time += realtime
                    cached_count += 1
                else:
                    executed_time += realtime
                    executed_count += 1
    
    except FileNotFoundError:
        print(f"Error: File '{trace_file}' not found")
        sys.exit(1)
    except Exception as e:
        print(f"Error reading trace file: {e}")
        sys.exit(1)
    
    total_time = executed_time + cached_time
    total_count = executed_count + cached_count
    
    print("=" * 60)
    print("WORKFLOW TIME ANALYSIS")
    print("=" * 60)
    print()
    print(f"📊 TASK SUMMARY:")
    print(f"  Tasks executed:    {executed_count:>6}")
    print(f"  Tasks cached:      {cached_count:>6}")
    print(f"  Total tasks:       {total_count:>6}")
    print()
    print(f"⏱️  TIME BREAKDOWN:")
    print(f"  Time executed:     {format_time(executed_time):>20}")
    print(f"  Time saved (cache): {format_time(cached_time):>20}")
    print()
    print("=" * 60)
    print(f"🕐 TOTAL TIME WITHOUT CACHE: {format_time(total_time)}")
    print("=" * 60)
    print()
    
    if total_time > 0:
        efficiency = (cached_time / total_time) * 100
        print(f"⚡ Cache efficiency: {efficiency:.1f}% time saved")
        
        if cached_count > 0:
            print(f"💾 Average cached task time: {format_time(cached_time / cached_count)}")
        if executed_count > 0:
            print(f"⚙️  Average executed task time: {format_time(executed_time / executed_count)}")
    
    print()

if __name__ == '__main__':
    if len(sys.argv) > 1:
        trace_file = sys.argv[1]
    else:
        # Try common trace file names
        import os
        possible_names = ['trace.txt', 'trace.tsv', 'execution_trace.txt']
        trace_file = None
        for name in possible_names:
            if os.path.exists(name):
                trace_file = name
                break
        
        if not trace_file:
            print("Usage: python3 analyze_trace.py [trace_file]")
            print("\nNo trace file specified and no default trace file found.")
            print(f"Looked for: {', '.join(possible_names)}")
            sys.exit(1)
        
        print(f"Using trace file: {trace_file}\n")
    
    analyze_trace(trace_file)
# Clean script for the React Rod Simulation component
Write-Host "🧹 Cleaning project artifacts..." -ForegroundColor Yellow

# Function to remove directory with retry logic
function Remove-DirectoryWithRetry {
    param(
        [string]$Path,
        [int]$MaxRetries = 3
    )
    
    if (-not (Test-Path $Path)) {
        Write-Host "✅ $Path does not exist, skipping" -ForegroundColor Green
        return
    }
    
    for ($i = 1; $i -le $MaxRetries; $i++) {
        try {
            Write-Host "Attempting to remove $Path (attempt $i/$MaxRetries)..." -ForegroundColor Cyan
            Remove-Item -Path $Path -Recurse -Force -ErrorAction Stop
            Write-Host "✅ Successfully removed $Path" -ForegroundColor Green
            return
        }
        catch {
            Write-Host "⚠️ Failed to remove $Path on attempt $i" -ForegroundColor Yellow
            if ($i -eq $MaxRetries) {
                Write-Host "❌ Could not remove $Path after $MaxRetries attempts. Some files may be locked." -ForegroundColor Red
                Write-Host "Try closing VS Code and any other applications that might be using these files." -ForegroundColor Red
            } else {
                Start-Sleep -Seconds 2
            }
        }
    }
}

# Remove node_modules
Remove-DirectoryWithRetry -Path "node_modules"

# Remove lockfiles
@("pnpm-lock.yaml", "package-lock.json", "yarn.lock") | ForEach-Object {
    if (Test-Path $_) {
        try {
            Remove-Item $_ -Force
            Write-Host "✅ Removed $_" -ForegroundColor Green
        }
        catch {
            Write-Host "⚠️ Could not remove $_" -ForegroundColor Yellow
        }
    }
}

# Remove dist directory
Remove-DirectoryWithRetry -Path "dist"

Write-Host "🎉 Clean operation completed!" -ForegroundColor Green

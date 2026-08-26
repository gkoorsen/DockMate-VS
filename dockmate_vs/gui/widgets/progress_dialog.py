"""
Progress dialog for campaign execution.

Shows progress bar, log output, and elapsed time during campaign execution.
"""

import tkinter as tk
from tkinter import ttk
import time
from datetime import datetime


class ProgressDialog(tk.Toplevel):
    """Modal progress dialog with log output."""

    def __init__(self, parent, total_ligands: int = 0):
        """
        Initialize progress dialog.

        Args:
            parent: Parent window
            total_ligands: Total number of ligands to dock
        """
        super().__init__(parent)

        self.title("DockMate-VS Campaign Running...")
        self.geometry("700x500")
        self.transient(parent)
        self.grab_set()  # Modal

        # Center on parent
        self.update_idletasks()
        x = parent.winfo_x() + (parent.winfo_width() // 2) - (700 // 2)
        y = parent.winfo_y() + (parent.winfo_height() // 2) - (500 // 2)
        self.geometry(f"700x500+{x}+{y}")

        self.total_ligands = total_ligands
        self.current_ligand = 0
        self.start_time = time.time()
        self.cancelled = False

        # Configure grid
        self.grid_rowconfigure(2, weight=1)
        self.grid_columnconfigure(0, weight=1)

        # Status label
        self.status_label = tk.Label(
            self,
            text="Initializing campaign...",
            font=('Arial', 11)
        )
        self.status_label.grid(row=0, column=0, pady=(20, 10), padx=20, sticky='w')

        # Progress bar
        progress_frame = tk.Frame(self)
        progress_frame.grid(row=1, column=0, sticky='ew', padx=20, pady=(0, 10))
        progress_frame.grid_columnconfigure(0, weight=1)

        self.progress = ttk.Progressbar(
            progress_frame,
            mode='determinate',
            maximum=100
        )
        self.progress.grid(row=0, column=0, sticky='ew')

        self.progress_label = tk.Label(
            progress_frame,
            text="0%",
            font=('Arial', 9),
            width=10
        )
        self.progress_label.grid(row=0, column=1, padx=(10, 0))

        # Log output frame
        log_frame = tk.LabelFrame(self, text="Log Output", padx=5, pady=5)
        log_frame.grid(row=2, column=0, sticky='nsew', padx=20, pady=(0, 10))
        log_frame.grid_rowconfigure(0, weight=1)
        log_frame.grid_columnconfigure(0, weight=1)

        self.log_text = tk.Text(
            log_frame,
            height=15,
            state='disabled',
            font=('Courier', 9),
            wrap='word'
        )
        self.log_text.grid(row=0, column=0, sticky='nsew')

        # Scrollbar for log
        log_scrollbar = ttk.Scrollbar(log_frame, orient='vertical', command=self.log_text.yview)
        log_scrollbar.grid(row=0, column=1, sticky='ns')
        self.log_text.configure(yscrollcommand=log_scrollbar.set)

        # Time labels
        time_frame = tk.Frame(self)
        time_frame.grid(row=3, column=0, sticky='ew', padx=20, pady=(0, 10))

        self.elapsed_label = tk.Label(
            time_frame,
            text="Elapsed: 0:00",
            font=('Arial', 9)
        )
        self.elapsed_label.pack(side='left')

        self.remaining_label = tk.Label(
            time_frame,
            text="Remaining: --:--",
            font=('Arial', 9)
        )
        self.remaining_label.pack(side='right')

        # Cancel button
        self.cancel_button = tk.Button(
            self,
            text="Cancel Campaign",
            command=self.cancel,
            width=15,
            bg='#e74c3c',
            fg='white'
        )
        self.cancel_button.grid(row=4, column=0, pady=(0, 20))

        # Start time update
        self._update_time()

    def update_progress(self, current: int, total: int, message: str):
        """
        Update progress bar and status.

        Args:
            current: Current progress value
            total: Total progress value
            message: Status message
        """
        self.current_ligand = current

        # Update progress bar
        if total > 0:
            percent = (current / total) * 100
            self.progress['value'] = percent
            self.progress_label.config(text=f"{percent:.0f}%")
        else:
            self.progress['mode'] = 'indeterminate'
            self.progress.start()

        # Update status
        self.status_label.config(text=message)

        # Log message
        self.log(f"[{current}/{total}] {message}")

        # Update display
        self.update_idletasks()

    def log(self, message: str):
        """
        Append message to log.

        Args:
            message: Log message
        """
        timestamp = datetime.now().strftime("%H:%M:%S")
        self.log_text.config(state='normal')
        self.log_text.insert('end', f"[{timestamp}] {message}\n")
        self.log_text.see('end')
        self.log_text.config(state='disabled')

    def _update_time(self):
        """Update elapsed and remaining time labels."""
        if self.cancelled:
            return

        # Calculate elapsed time
        elapsed_seconds = time.time() - self.start_time
        elapsed_minutes = int(elapsed_seconds // 60)
        elapsed_secs = int(elapsed_seconds % 60)
        self.elapsed_label.config(text=f"Elapsed: {elapsed_minutes}:{elapsed_secs:02d}")

        # Calculate remaining time (rough estimate)
        if self.current_ligand > 0 and self.total_ligands > 0:
            time_per_ligand = elapsed_seconds / self.current_ligand
            remaining_ligands = self.total_ligands - self.current_ligand
            remaining_seconds = time_per_ligand * remaining_ligands
            remaining_minutes = int(remaining_seconds // 60)
            remaining_secs = int(remaining_seconds % 60)
            self.remaining_label.config(text=f"Remaining: ~{remaining_minutes}:{remaining_secs:02d}")
        else:
            self.remaining_label.config(text="Remaining: --:--")

        # Schedule next update
        self.after(1000, self._update_time)

    def cancel(self):
        """Cancel campaign."""
        result = tk.messagebox.askyesno(
            "Cancel Campaign",
            "Are you sure you want to cancel the campaign?\n\n"
            "This will stop the docking process."
        )
        if result:
            self.cancelled = True
            self.log("Campaign cancelled by user")
            self.destroy()

    def complete(self):
        """Mark campaign as complete."""
        self.status_label.config(text="Campaign completed successfully!")
        self.progress['value'] = 100
        self.progress_label.config(text="100%")
        self.cancel_button.config(
            text="Close",
            bg='#27ae60',
            command=self.destroy
        )
        self.log("Campaign completed!")

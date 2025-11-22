"""plotting.py

This module defines functions for plotting wake behavior and body force coeffs

Functions:
- gridtransform(...)
- update_colorbars(...)
- progress_callback(...)
- animation_saving(...)
"""

import numpy as np
import matplotlib
import os


def gridtransform(gs, scaling=1.0, x_offset=0.0, y_offset=0.0, elements=None):
    """Applies postprocess scaling for the plot based on a scaling factor and offset
    """
    # Collect element positional values
    x0_list, x1_list, y0_list, y1_list = [], [], [], []
    for ax in elements:
        pos = ax.get_position()
        x0_list.append(pos.x0)
        x1_list.append(pos.x1)
        y0_list.append(pos.y0)
        y1_list.append(pos.y1)
    
    # Original Bounding Box
    bbox = np.array([min(x0_list), max(x1_list), min(y0_list), max(y1_list)], dtype=float)
    
    # Center of the bbox
    center_x = (bbox[0] + bbox[1]) / 2
    center_y = (bbox[2] + bbox[3]) / 2
    
    # Scale about center and apply offset
    new_left   = (bbox[0] - center_x) * scaling + center_x + x_offset
    new_right  = (bbox[1] - center_x) * scaling + center_x + x_offset
    new_bottom = (bbox[2] - center_y) * scaling + center_y + y_offset
    new_top    = (bbox[3] - center_y) * scaling + center_y + y_offset
    
    # Apply to GridSpec
    gs.update(left=new_left, right=new_right, bottom=new_bottom, top=new_top)


def update_colorbars(fig, right_ax, under_ax, boundplot, wakeplot,
                     norm_bound, norm_wake, mode='mixed'):
    """
    Adjust colorbar visibility, positions, and normalization based on selected mode.

    Modes:
    'boundscaled'   :   color scale for vortex strength is scaled against bound vortex strength
    'wakescaled'    :   color scale for vortex strength is scaled against wake vortex strength
    'mixed'         :   scales wake and bound vortices independently on separate color bars (preferred)
    """
    
    # Default: hide both colorbars
    right_ax.set_visible(False)
    under_ax.set_visible(False)
    cbar_right, cbar_under = None, None
    
    if mode == 'boundscaled':
        right_ax.set_visible(True)
        boundplot.set_norm(norm_bound)
        wakeplot.set_norm(norm_bound)
        # Initialize colorbar
        cbar_right = fig.colorbar(boundplot, cax=right_ax)
        cbar_right.set_label("Vortex Strength Γ")
        
    elif mode == 'wakescaled':
        right_ax.set_visible(True)
        boundplot.set_norm(norm_wake)
        wakeplot.set_norm(norm_wake)
        # Initialize colorbar
        cbar_right = fig.colorbar(wakeplot, cax=right_ax)
        cbar_right.set_label("Vortex Strength Γ")
        
    elif mode == 'mixed':
        # Show both colorbars
        right_ax.set_visible(True)
        under_ax.set_visible(True)
        # Normalize both scatters
        boundplot.set_norm(norm_bound)
        wakeplot.set_norm(norm_wake)
        # Initialize colorbars
        cbar_right = fig.colorbar(boundplot, cax=right_ax)
        cbar_right.set_label("Bound Vortex Γ")
        cbar_under = fig.colorbar(wakeplot, cax=under_ax, orientation='horizontal')
        cbar_under.ax.xaxis.set_ticks_position('bottom')
        cbar_under.ax.xaxis.set_label_position('bottom')
        cbar_under.set_label("Wake Vortex Γ", rotation=0, labelpad=8)
        
    else:
        raise ValueError(f"Invalid mode '{mode}'. Must be 'boundscaled', 'wakescaled', or 'mixed'.")
    
    if cbar_right:
        cbar_right.update_normal(boundplot if mode != 'wakescaled' else wakeplot)
    if cbar_under:
        cbar_under.update_normal(wakeplot)



def progress_callback(current_frame, total_frames):
    """Display progress of animation saving in the console."""
    print("\r", end="", flush=True)
    if current_frame == (total_frames - 1):
        print("Encoding file...          ", end="\n", flush=True)
    else:
        print(f"Saving frame {current_frame}/{total_frames}", end="", flush=True)


def animation_saving(ani, format, file_name, path=".", fps=20, dpi='figure', bitrate=-1, embed_limit=50, progress=progress_callback):
    """Save the animated plot as varous formats.

    format : str
        "html5video" : This saves the animation as an h264 video, encoded in base64 directly into the HTML5 video tag. Uses the FFMpegWriter
        "jshtml" : An HTML representation of the animation embedded as a js object. Uses the HTMLWriter
        "gif" : Saves the animation as a gif using the PillowWriter
        "mp4" : Saves the animation as a mp4 using the FFMpegWriter
        "html" : Saves the animation as a standalone html document with the individual frames stored in a supporting folder
    """

    # Makes the path if it doesnt exist already
    os.makedirs(path, exist_ok=True)

    # Construct full file path 
    full_path = lambda ext: os.path.join(path, f"{file_name}.{ext}")

    print("Saving Animation...")
    
    match format:
        case "html5video":
            try:
                import imageio_ffmpeg
            except ImportError:
                print("FFmpeg not installed. FFmpeg is required for 'html5video' and 'mp4' saving methods, try installing package 'imageio_ffmpeg'.")
                return
            else:
                matplotlib.rcParams["animation.ffmpeg_path"] = imageio_ffmpeg.get_ffmpeg_exe()

            print("Standby...")
            matplotlib.rcParams['animation.bitrate'] = bitrate
            html_str = ani.to_html5_video(embed_limit=embed_limit)
            with open(full_path("html"), "w", encoding="utf-8") as f:
                f.write(html_str) 

        case "jshtml":
            print("Standby...")
            matplotlib.rcParams['animation.embed_limit'] = embed_limit  # in MB
            html_str = ani.to_jshtml(fps=fps, embed_frames=True, default_mode='loop')
            with open(full_path("html"), "w", encoding="utf-8") as f:
                f.write(html_str)

        case "gif":
            ani.save(filename=full_path("gif"), writer="pillow", fps=fps, dpi=dpi, bitrate=bitrate, progress_callback=progress) 
        
        case "mp4":
            try:
                import imageio_ffmpeg
            except ImportError:
                print("FFmpeg not installed. FFmpeg is required for 'html5video' and 'mp4' saving methods, try installing package 'imageio_ffmpeg'.")
                return
            else:
                matplotlib.rcParams["animation.ffmpeg_path"] = imageio_ffmpeg.get_ffmpeg_exe()

            ani.save(filename=full_path("mp4"), writer="ffmpeg", fps=fps, dpi=dpi, bitrate=bitrate, progress_callback=progress) 
        
        case "html":
            ani.save(filename=full_path("html"), writer="html", fps=fps, dpi=dpi, bitrate=bitrate, progress_callback=progress) 
    
    print(f"Animation saved to: {os.path.abspath(path)}")
    
        


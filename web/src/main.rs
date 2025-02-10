use bsxplorer_ci::convert::{ConvertReportType, IpcCompression};
use gloo_timers::callback::Interval;
use stylist::style;
use stylist::yew::styled_component;
use wasm_bindgen::prelude::wasm_bindgen;
use wasm_bindgen::JsValue;
use wasm_bindgen_futures::JsFuture;
use web_sys::{FileSystemDirectoryHandle, HtmlInputElement};
use yew;
use yew::events::SubmitEvent;
use yew::prelude::*;
use yew_router::prelude::*;

#[derive(Clone, Routable, PartialEq)]
enum Route {
    #[at("/")]
    Home,
    #[at("/about")]
    About,
    #[at("/dmr")]
    Dmr,
    #[at("/convert")]
    Convert,
    #[at("/stats")]
    Stats,
    #[not_found]
    #[at("/404")]
    NotFound,
}

#[function_component(Home)]
fn home() -> Html {
    html! {
        <div>
            <h1>{ "Home Page" }</h1>
            <p>{ "Welcome to the Home page." }</p>
        </div>
    }
}

#[function_component(About)]
fn about() -> Html {
    html! {
        <div>
            <h1>{ "About Page" }</h1>
            <p>{ "Learn more About us on this page." }</p>
        </div>
    }
}

#[function_component(Stats)]
fn stats() -> Html {
    html! {
        <div>
            <h1>{ "About Page" }</h1>
            <p>{ "Learn more About us on this page." }</p>
        </div>
    }
}

#[wasm_bindgen]
extern "C" {
    #[wasm_bindgen(js_name = "showDirectoryPicker")]
    fn show_directory_picker() -> JsValue;
}

#[derive(Clone, Debug, PartialEq)]
pub struct Config {
    pub input_files: Vec<String>,
    pub output_dir: String,
    pub from_type: ConvertReportType,
    pub into_type: ConvertReportType,
    pub ipc_compression: Option<IpcCompression>,
    pub threads: u32,
    pub low_memory: bool,
    pub chunk_size: u32,
    pub fasta_path: Option<String>,
    pub fai_path: Option<String>,
    pub batch_per_read: u32,
    pub batch_size: u32,
}

impl Default for Config {
    fn default() -> Self {
        Self {
            input_files: vec![],
            output_dir: "".into(),
            from_type: ConvertReportType::Bsx,
            into_type: ConvertReportType::Bsx,
            ipc_compression: Some(IpcCompression::None),
            // In a real app you might detect the number of system threads;
            // here we default to 4 for demonstration purposes.
            threads: 4,
            low_memory: false,
            chunk_size: 20000,
            fasta_path: None,
            fai_path: None,
            batch_per_read: 8,
            batch_size: 1000000,
        }
    }
}

#[derive(Clone, Debug, PartialEq)]
enum ViewState {
    Config,
    Confirm,
    Progress,
    Done,
}

#[styled_component(Convert)]
fn convert() -> Html {
    let config = use_state(|| Config::default());
    let view_state = use_state(|| ViewState::Config);
    let progress = use_state(|| 0u32);

    // --- File input for multiple input files ---
    let on_input_files_change = {
        let config = config.clone();
        Callback::from(move |e: Event| {
            let input: HtmlInputElement = e.target_unchecked_into();
            let files = input.files();
            let mut file_names = vec![];
            if let Some(files_list) = files {
                for i in 0..files_list.length() {
                    if let Some(file) = files_list.get(i) {
                        file_names.push(file.name());
                    }
                }
            }
            let mut new_config = (*config).clone();
            new_config.input_files = file_names;
            config.set(new_config);
        })
    };

    let output_dir = use_state(|| String::new());

    let on_pick_directory = {
        let output_dir = output_dir.clone();
        Callback::from(move |e: MouseEvent| {
            // Capture the click event
            e.prevent_default(); // Prevent form submission
            e.stop_propagation(); // Stop event bubbling

            let output_dir_clone = output_dir.clone();

            let promise = show_directory_picker();
            wasm_bindgen_futures::spawn_local(async move {
                if let Ok(js_value) = JsFuture::from(js_sys::Promise::from(promise)).await {
                    let dir_handle: FileSystemDirectoryHandle = js_value.into();
                    let dir_name = dir_handle.name();
                    output_dir_clone.set(dir_name);
                }
            });
        })
    };

    // --- Output directory selection (using webkitdirectory) ---
    let on_output_dir_change = {
        let config = config.clone();
        Callback::from(move |e: Event| {
            // For a directory selection input, we extract one file's path.
            let input: HtmlInputElement = e.target_unchecked_into();
            let files = input.files();
            let dir_name = if let Some(files_list) = files {
                if files_list.length() > 0 {
                    if let Some(file) = files_list.get(0) {
                        // Here we use the file name as a placeholder.
                        file.name()
                    } else {
                        "".into()
                    }
                } else {
                    "".into()
                }
            } else {
                "".into()
            };
            let mut new_config = (*config).clone();
            new_config.output_dir = dir_name;
            config.set(new_config);
        })
    };

    // --- Selection for from_type ---
    let on_from_type_change = {
        let config = config.clone();
        Callback::from(move |e: Event| {
            let select: web_sys::HtmlSelectElement = e.target_unchecked_into();
            let value = select.value();
            let from_type = match value.as_str() {
                "Bsx" => ConvertReportType::Bsx,
                "Bismark" => ConvertReportType::Bismark,
                "CgMap" => ConvertReportType::CgMap,
                "BedGraph" => ConvertReportType::BedGraph,
                "Coverage" => ConvertReportType::Coverage,
                _ => ConvertReportType::Bsx,
            };
            let mut new_config = (*config).clone();
            new_config.from_type = from_type;
            config.set(new_config);
        })
    };

    // --- Selection for into_type ---
    let on_into_type_change = {
        let config = config.clone();
        Callback::from(move |e: Event| {
            let select: web_sys::HtmlSelectElement = e.target_unchecked_into();
            let value = select.value();
            let into_type = match value.as_str() {
                "Bsx" => ConvertReportType::Bsx,
                "Bismark" => ConvertReportType::Bismark,
                "CgMap" => ConvertReportType::CgMap,
                "BedGraph" => ConvertReportType::BedGraph,
                "Coverage" => ConvertReportType::Coverage,
                _ => ConvertReportType::Bsx,
            };
            let mut new_config = (*config).clone();
            new_config.into_type = into_type;
            // If the selected into_type is not Bsx, clear ipc_compression.
            if let ConvertReportType::Bsx = new_config.into_type {
                if new_config.ipc_compression.is_none() {
                    new_config.ipc_compression = Some(IpcCompression::None);
                }
            } else {
                new_config.ipc_compression = None;
            }
            config.set(new_config);
        })
    };

    // --- Selection for ipc_compression (only visible if into_type is Bsx) ---
    let on_ipc_compression_change = {
        let config = config.clone();
        Callback::from(move |e: Event| {
            let select: web_sys::HtmlSelectElement = e.target_unchecked_into();
            let value = select.value();
            let ipc = match value.as_str() {
                "LZ4" => IpcCompression::LZ4,
                "ZSTD" => IpcCompression::ZSTD,
                "None" => IpcCompression::None,
                _ => IpcCompression::None,
            };
            let mut new_config = (*config).clone();
            new_config.ipc_compression = Some(ipc);
            config.set(new_config);
        })
    };

    // --- Threads number input ---
    let on_threads_change = {
        let config = config.clone();
        Callback::from(move |e: Event| {
            let input: HtmlInputElement = e.target_unchecked_into();
            if let Ok(num) = input.value().parse::<u32>() {
                let mut new_config = (*config).clone();
                new_config.threads = num;
                config.set(new_config);
            }
        })
    };

    // --- Low memory checkbox ---
    let on_low_memory_change = {
        let config = config.clone();
        Callback::from(move |e: Event| {
            let input: HtmlInputElement = e.target_unchecked_into();
            let checked = input.checked();
            let mut new_config = (*config).clone();
            new_config.low_memory = checked;
            config.set(new_config);
        })
    };

    // --- Chunk size number input ---
    let on_chunk_size_change = {
        let config = config.clone();
        Callback::from(move |e: Event| {
            let input: HtmlInputElement = e.target_unchecked_into();
            if let Ok(num) = input.value().parse::<u32>() {
                let mut new_config = (*config).clone();
                new_config.chunk_size = num;
                config.set(new_config);
            }
        })
    };

    // --- FASTA path: required if from_type is BedGraph or Coverage ---
    let on_fasta_path_change = {
        let config = config.clone();
        Callback::from(move |e: Event| {
            let input: HtmlInputElement = e.target_unchecked_into();
            let files = input.files();
            let fasta = if let Some(files_list) = files {
                if files_list.length() > 0 {
                    if let Some(file) = files_list.get(0) {
                        Some(file.name())
                    } else {
                        None
                    }
                } else {
                    None
                }
            } else {
                None
            };
            let mut new_config = (*config).clone();
            new_config.fasta_path = fasta;
            config.set(new_config);
        })
    };

    // --- FAI path: required if from_type is BedGraph/Coverage or into_type is Bsx ---
    let on_fai_path_change = {
        let config = config.clone();
        Callback::from(move |e: Event| {
            let input: HtmlInputElement = e.target_unchecked_into();
            let files = input.files();
            let fai = if let Some(files_list) = files {
                if files_list.length() > 0 {
                    if let Some(file) = files_list.get(0) {
                        Some(file.name())
                    } else {
                        None
                    }
                } else {
                    None
                }
            } else {
                None
            };
            let mut new_config = (*config).clone();
            new_config.fai_path = fai;
            config.set(new_config);
        })
    };

    // --- Batch per read number input ---
    let on_batch_per_read_change = {
        let config = config.clone();
        Callback::from(move |e: Event| {
            let input: HtmlInputElement = e.target_unchecked_into();
            if let Ok(num) = input.value().parse::<u32>() {
                let mut new_config = (*config).clone();
                new_config.batch_per_read = num;
                config.set(new_config);
            }
        })
    };

    // --- Batch size number input ---
    let on_batch_size_change = {
        let config = config.clone();
        Callback::from(move |e: Event| {
            let input: HtmlInputElement = e.target_unchecked_into();
            if let Ok(num) = input.value().parse::<u32>() {
                let mut new_config = (*config).clone();
                new_config.batch_size = num;
                config.set(new_config);
            }
        })
    };

    // ----------------------------
    // on_submit: Switch from config to confirm view.
    let on_submit = {
        let view_state = view_state.clone();
        Callback::from(move |e: SubmitEvent| {
            e.prevent_default();
            // In a complete app, you might validate & update `config` here.
            view_state.set(ViewState::Confirm);
        })
    };

    // ----------------------------
    // on_confirm: Start a dummy function to simulate progress.
    let on_confirm = {
        let view_state = view_state.clone();
        let progress = progress.clone();
        Callback::from(move |_| {
            view_state.set(ViewState::Progress);
            // Start an interval timer to update progress.
            let progress_clone = progress.clone();
            let view_state = view_state.clone();
            let interval = Interval::new(100, move || {
                let current = *progress_clone;
                if current < 100 {
                    progress_clone.set(current + 5);
                } else {
                    view_state.set(ViewState::Done);
                }
            });
            // Keep the interval alive.
            interval.forget();
        })
    };

    // ----------------------------
    // on_back: Return to the config screen.
    let on_back = {
        let view_state = view_state.clone();
        Callback::from(move |_| {
            view_state.set(ViewState::Config);
        })
    };
    let is_confirm_disabled = config.input_files.is_empty() || config.output_dir.is_empty();

    let input_files_class = if config.input_files.is_empty() {
        "error"
    } else {
        ""
    };
    let output_dir_class = if config.output_dir.is_empty() {
        "error"
    } else {
        ""
    };

    // --- Styling with stylist ---
    let stylesheet = stylist::style!(
        r#"
        .container {
            display: flex;
            flex-direction: column;
            gap: 1rem;
            max-width: 600px;
            margin: 2rem auto;
            padding: 1rem;
            background-color: #f9f9f9;
            border-radius: 8px;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
        }
        .field {
            display: flex;
            flex-direction: column;
        }
        .field label {
            margin-bottom: 0.5rem;
            font-weight: bold;
        }
        .field input, .field select {
            padding: 0.5rem;
            font-size: 1rem;
            border: 1px solid #ccc;
            border-radius: 4px;
        }
        .checkbox-field {
            flex-direction: row;
            align-items: center;
        }
        .checkbox-field input {
            margin-right: 0.5rem;
        }
        button {
            padding: 0.75rem;
            font-size: 1rem;
            background-color: #007BFF;
            color: white;
            border: none;
            border-radius: 4px;
            cursor: pointer;
        }
        button:hover {
            background-color: #0056b3;
        }
        submit::disabled {
            background-color: #a5d6a7;
            cursor: not-allowed;
        }
        .deny {
            background-color: #dc3545;
            color: white;
        }
        .error {
            border: 2px solid red;
            background-color: #ffe6e6;
        }
        .button-container {
            display: flex;
            justify-content: space-between;
            margin-top: 1rem;
        }
        .confirm {
            background-color: #28a745; /* green */
            color: white;
        }
        .deny {
            background-color: #dc3545; /* red */
            color: white;
        }
        "#
    )
    .unwrap();

    html! {
        <div class={stylesheet}>
            {
                match *view_state {
                    ViewState::Config => html! {
                        <form class="container" onsubmit={on_submit}>
                        // --- Input files (multiple) ---
                        <div class="field">
                            <label for="input_files">{"Input Files"}</label>
                            <input id="input_files" type="file" multiple=true onchange={on_input_files_change} class={input_files_class}/>
                        </div>

                        // --- Output directory selection ---
                        <div class="field">
                            <label for="output_dir">{ "Output Directory:" }</label>
                            <button onclick={on_pick_directory} class="select-directory">{ "Select Directory" }</button>
                            <p>{ format!("Selected: {}", *output_dir) }</p>
                        </div>

                        // --- From Type selection ---
                        <div class="field">
                            <label for="from_type">{"From Type"}</label>
                            <select id="from_type" onchange={on_from_type_change}>
                                <option value="Bsx" selected={config.from_type == ConvertReportType::Bsx}>{"Bsx"}</option>
                                <option value="Bismark" selected={config.from_type == ConvertReportType::Bismark}>{"Bismark"}</option>
                                <option value="CgMap" selected={config.from_type == ConvertReportType::CgMap}>{"CgMap"}</option>
                                <option value="BedGraph" selected={config.from_type == ConvertReportType::BedGraph}>{"BedGraph"}</option>
                                <option value="Coverage" selected={config.from_type == ConvertReportType::Coverage}>{"Coverage"}</option>
                            </select>
                        </div>

                        // --- Into Type selection ---
                        <div class="field">
                            <label for="into_type">{"Into Type"}</label>
                            <select id="into_type" onchange={on_into_type_change}>
                                <option value="Bsx" selected={config.into_type == ConvertReportType::Bsx}>{"Bsx"}</option>
                                <option value="Bismark" selected={config.into_type == ConvertReportType::Bismark}>{"Bismark"}</option>
                                <option value="CgMap" selected={config.into_type == ConvertReportType::CgMap}>{"CgMap"}</option>
                                <option value="BedGraph" selected={config.into_type == ConvertReportType::BedGraph}>{"BedGraph"}</option>
                                <option value="Coverage" selected={config.into_type == ConvertReportType::Coverage}>{"Coverage"}</option>
                            </select>
                        </div>

                        // --- IPC Compression (only shown if into_type == Bsx) ---
                        { if let ConvertReportType::Bsx = config.into_type {
                            html! {
                                <div class="field">
                                    <label for="ipc_compression">{"IPC Compression"}</label>
                                    <select id="ipc_compression" onchange={on_ipc_compression_change}>
                                        <option value="LZ4">{"LZ4"}</option>
                                        <option value="ZSTD">{"ZSTD"}</option>
                                        <option value="None" selected={matches!(config.ipc_compression, Some(IpcCompression::None))}>{"None"}</option>
                                    </select>
                                </div>
                            }
                        } else {
                            html! {}
                        }}

                        // --- Threads input ---
                        <div class="field">
                            <label for="threads">{"Threads"}</label>
                            <input
                                id="threads"
                                type="number"
                                min="1"
                                max="16"
                                value={config.threads.to_string()}
                                onchange={on_threads_change}
                            />
                        </div>

                        // --- Low Memory checkbox ---
                        <div class="field checkbox-field">
                            <input
                                id="low_memory"
                                type="checkbox"
                                checked={config.low_memory}
                                onchange={on_low_memory_change}
                            />
                            <label for="low_memory">{"Low Memory"}</label>
                        </div>

                        // --- Chunk Size input ---
                        <div class="field">
                            <label for="chunk_size">{"Chunk Size"}</label>
                            <input
                                id="chunk_size"
                                type="number"
                                value={config.chunk_size.to_string()}
                                onchange={on_chunk_size_change}
                            />
                        </div>

                        // --- FASTA Path: required if from_type is BedGraph or Coverage ---
                        { if matches!(config.from_type, ConvertReportType::BedGraph | ConvertReportType::Coverage) {
                            html! {
                                <div class="field">
                                    <label for="fasta_path">{"FASTA Path (required)"}</label>
                                    <input id="fasta_path" type="file" onchange={on_fasta_path_change} />
                                </div>
                            }
                        } else {
                            html! {
                                <div class="field">
                                    <label for="fasta_path">{"FASTA Path (optional)"}</label>
                                    <input id="fasta_path" type="file" onchange={on_fasta_path_change} />
                                </div>
                            }
                        }}

                        // --- FAI Path: required if from_type is BedGraph/Coverage or into_type is Bsx ---
                        { if matches!(config.from_type, ConvertReportType::BedGraph | ConvertReportType::Coverage) || matches!(config.into_type, ConvertReportType::Bsx) {
                            html! {
                                <div class="field">
                                    <label for="fai_path">{"FAI Path (required)"}</label>
                                    <input id="fai_path" type="file" onchange={on_fai_path_change} />
                                </div>
                            }
                        } else {
                            html! {
                                <div class="field">
                                    <label for="fai_path">{"FAI Path (optional)"}</label>
                                    <input id="fai_path" type="file" onchange={on_fai_path_change} />
                                </div>
                            }
                        }}

                        // --- Batch Per Read input ---
                        <div class="field">
                            <label for="batch_per_read">{"Batch Per Read"}</label>
                            <input
                                id="batch_per_read"
                                type="number"
                                value={config.batch_per_read.to_string()}
                                onchange={on_batch_per_read_change}
                            />
                        </div>

                        // --- Batch Size input ---
                        <div class="field">
                            <label for="batch_size">{"Batch Size"}</label>
                            <input
                                id="batch_size"
                                type="number"
                                value={config.batch_size.to_string()}
                                onchange={on_batch_size_change}
                            />
                        </div>

                        <button disabled={is_confirm_disabled} type="submit" class="submit">{"Submit"}</button>
                    </form>
                },

                ViewState::Confirm => {
                    html! {
                        <div class="container">
                            <h2>{ "Confirm Your Configuration" }</h2>
                            <p>{ format!("Input Files: {:?}", config.input_files) }</p>
                            <p>{ format!("Output Directory: {}", config.output_dir) }</p>
                            <p>{ format!("From Type: {:?}", config.from_type) }</p>
                            <p>{ format!("Into Type: {:?}", config.into_type) }</p>
                            // ... Display other configuration parameters ...
                            <div class="button-container">
                                <button class="deny" onclick={on_back.clone()}>{ "Deny" }</button>
                                <button class="confirm" onclick={on_confirm.clone()}>{ "Confirm" }</button>
                            </div>
                        </div>
                    }
                },

                ViewState::Progress => {
                    html! {
                        <div class="container">
                            <h2>{ "Processing..." }</h2>
                            <progress max="100" value={progress.to_string()}>{ format!("{}%", *progress) }</progress>
                        </div>
                    }
                },
                ViewState::Done => {
                    html! {
                        <div class="container">
                            <h2>{ "Process Complete" }</h2>
                            <button onclick={on_back.clone()}>{ "Back to Config" }</button>
                        </div>
                    }
                },
            }
        }
        </div>
    }
}

#[function_component(Dmr)]
fn dmr() -> Html {
    html! {
        <div>
            <h1>{ "About Page" }</h1>
            <p>{ "Learn more About us on this page." }</p>
        </div>
    }
}

#[function_component(NavMenu)]
fn nav_menu() -> Html {
    let is_menu_open = use_state(|| false);
    let is_settings_open = use_state(|| false);

    let toggle_menu = {
        let is_menu_open = is_menu_open.clone();
        Callback::from(move |_| is_menu_open.set(!*is_menu_open))
    };

    let toggle_settings = {
        let is_settings_open = is_settings_open.clone();
        Callback::from(move |_| is_settings_open.set(!*is_settings_open))
    };

    html! {
        <nav class="navbar">
            <div class="menu-toggle" onclick={toggle_menu.clone()}>
                <div class="bar"></div>
                <div class="bar"></div>
                <div class="bar"></div>
            </div>

            // Navigation Links
            <ul class={if *is_menu_open { "nav-links open" } else { "nav-links" }}>
                <li><Link<Route> to={Route::Home}>{ "Home" }</Link<Route>></li>
                <li><Link<Route> to={Route::About}>{ "About" }</Link<Route>></li>
                <li><Link<Route> to={Route::Dmr}>{ "DMR" }</Link<Route>></li>
                <li><Link<Route> to={Route::Convert}>{ "Convert" }</Link<Route>></li>
                <li><Link<Route> to={Route::Stats}>{ "Stats" }</Link<Route>></li>
            </ul>

            // Settings Button (Top-Right Corner)
            <div class="settings-container">
                <button class="settings-btn" onclick={toggle_settings.clone()}>{ "⚙️" }</button>
                {
                    if *is_settings_open {
                        html! {
                            <div class="settings-menu">
                                <ul>
                                    <li>{ "Profile" }</li>
                                </ul>
                            </div>
                        }
                    } else {
                        html! {}
                    }
                }
            </div>
        </nav>
    }
}

#[function_component(App)]
fn app() -> Html {
    html! {
        <BrowserRouter>
            <NavMenu />
            <Switch<Route> render={move |routes| match routes {
                Route::Home => html! { <Home /> },
                Route::About => html! { <About /> },
                Route::Dmr => html! { <Dmr /> },
                Route::Convert => html! { <Convert /> },
                Route::Stats => html! { <Stats /> },
                Route::NotFound => html! { <h1>{ "404 - Page Not Found" }</h1> },
            } } />
        </BrowserRouter>
    }
}

fn main() {
    yew::Renderer::<App>::new().render();
}
